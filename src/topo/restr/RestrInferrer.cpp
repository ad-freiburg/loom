// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include <unordered_map>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/restr/RestrGraph.h"
#include "topo/restr/RestrInferrer.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Dijkstra.h"

using topo::restr::RestrInferrer;

// _____________________________________________________________________________
RestrInferrer::RestrInferrer(const TopoConfig* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void RestrInferrer::init(TransitGraph* g) {
  for (auto nd : *g->getNds()) {
    _nMap[nd] = _g.addNd();
  }

  for (auto nd : *g->getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      _eMap[edg] = {_g.addEdg(_nMap[edg->getFrom()], _nMap[edg->getTo()],
                              edg->pl().getPolyline()),
                    _g.addEdg(_nMap[edg->getTo()], _nMap[edg->getFrom()],
                              edg->pl().getPolyline().getReversed())};
      for (auto r : edg->pl().getRoutes()) {
        if (r.direction == 0) {
          _eMap[edg][0]->pl().routes.insert(r.route);
          _eMap[edg][1]->pl().routes.insert(r.route);
        }

        if (r.direction == edg->getTo()) {
          _eMap[edg][0]->pl().routes.insert(r.route);
        }

        if (r.direction == edg->getFrom()) {
          _eMap[edg][1]->pl().routes.insert(r.route);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void RestrInferrer::infer(TransitGraph* g, const OrigEdgs& origEdgs) {
  std::cerr << "Inserting handles..." << std::endl;

  insertHandles(g, origEdgs);

  util::geo::output::GeoGraphJsonOutput out;
  std::ofstream of;
  of.open("restr.graph");
  out.print(_g, of);
  of.flush();
  //
  std::cerr << "Inferring restrictons..." << std::endl;

  for (auto nd : *g->getNds()) {
    for (auto edg1 : nd->getAdjList()) {
      // check every other edge
      for (auto edg2 : nd->getAdjList()) {
        if (edg1 == edg2) continue;

        for (auto ro1 : edg1->pl().getRoutes()) {
          if (!edg2->pl().hasRoute(ro1.route)) continue;

          const auto& ro2 = edg2->pl().getRouteOcc(ro1.route);

          if (ro1.direction != 0 && ro2.direction != 0 &&
              ro1.direction == ro2.direction)
            continue;

          if (ro1.direction != 0 && ro2.direction != 0 &&
              edg1->getOtherNd(ro1.direction) ==
                  edg2->getOtherNd(ro2.direction)) {
            continue;
          }

          if (!check(ro1.route, edg1, edg2, origEdgs) &&
              !check(ro1.route, edg2, edg1, origEdgs)) {
            std::cerr << " ==> INFER" << std::endl;
            nd->pl().addConnExc(ro1.route, edg1, edg2);
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void RestrInferrer::insertHandles(TransitGraph* g, const OrigEdgs& origEdgs) {
  std::map<RestrEdge*, std::vector<std::pair<RestrNode*, double>>> handles;

  for (auto nd : *g->getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      insertHandles(edg, origEdgs, &handles);
    }
  }

  struct {
    bool operator()(const std::pair<RestrNode*, double>& a,
                    const std::pair<RestrNode*, double>& b) const {
      return a.second < b.second;
    }
  } cmp;

  for (auto edgHndl : handles) {
    std::sort(edgHndl.second.begin(), edgHndl.second.end(), cmp);

    auto lastNd = edgHndl.first->getFrom();
    double lastPos = 0;
    for (auto hndl : edgHndl.second) {
      std::cerr << hndl.second << std::endl;
      auto e =
          _g.addEdg(lastNd, hndl.first,
                    edgHndl.first->pl().geom.getSegment(lastPos, hndl.second));
      e->pl().routes = edgHndl.first->pl().routes;
      lastNd = hndl.first;
      lastPos = hndl.second;
    }
    auto e = _g.addEdg(lastNd, edgHndl.first->getTo(),
                       edgHndl.first->pl().geom.getSegment(lastPos, 1));
    e->pl().routes = edgHndl.first->pl().routes;

    _g.delEdg(edgHndl.first->getFrom(), edgHndl.first->getTo());
  }
}

// _____________________________________________________________________________
void RestrInferrer::insertHandles(
    const TransitEdge* e, const OrigEdgs& origEdgs,
    std::map<RestrEdge*, std::vector<std::pair<RestrNode*, double>>>* handles) {
  // TODO: use max_aggr_dist here
  double MAX_DIST = 60;

  auto handlePointA = e->pl().getPolyline().getPointAt(0.33).p;
  auto handlePointB = e->pl().getPolyline().getPointAt(0.66).p;

  auto handleLineA = e->pl().getPolyline().getOrthoLineAtDist(0.33 * e->pl().getPolyline().getLength(), e->pl().getRoutes().size() * MAX_DIST).getLine();
  auto handleLineB = e->pl().getPolyline().getOrthoLineAtDist(0.66 * e->pl().getPolyline().getLength(), e->pl().getRoutes().size() * MAX_DIST).getLine();

  for (auto edg : origEdgs.find(e)->second) {
    auto origFr = const_cast<TransitEdge*>(edg);

    const auto& edgs = _eMap.find(origFr)->second;


    for (auto restrE : edgs) {
      if (util::geo::intersects(handleLineA, restrE->pl().geom.getLine())) {
        auto projA = restrE->pl().geom.projectOn(handlePointA);
        auto handleNdA = _g.addNd();

        _handlesA[e].insert(handleNdA);

        (*handles)[restrE].push_back({handleNdA, projA.totalPos});
      }

      if (util::geo::intersects(handleLineB, restrE->pl().geom.getLine())) {
        auto projB = restrE->pl().geom.projectOn(handlePointB);
        auto handleNdB = _g.addNd();

        _handlesB[e].insert(handleNdB);

        (*handles)[restrE].push_back({handleNdB, projB.totalPos});
      }
    }
  }
}

// _____________________________________________________________________________
bool RestrInferrer::check(const Route* r, const TransitEdge* edg1,
                          const TransitEdge* edg2,
                          const OrigEdgs& origEdgs) const {
  std::set<RestrEdge*> from;
  std::set<RestrEdge*> to;

  std::map<const RestrEdge*, double> sourcePos;
  std::map<const RestrEdge*, double> targetPos;

  auto shrdNd = shared::transitgraph::TransitGraph::sharedNode(edg1, edg2);

  double curDist = edg1->pl().getPolyline().getLength() * 0.33 +
                   edg2->pl().getPolyline().getLength() * 0.33;

  if (shrdNd == edg1->getFrom()) {
    if (_handlesA.count(edg1)) {
      for (auto nd : _handlesA.find(edg1)->second) {
        from.insert(nd->getAdjListIn().begin(), nd->getAdjListIn().end());
      }
    }
  } else {
    if (_handlesB.count(edg1)) {
      for (auto nd : _handlesB.find(edg1)->second) {
        from.insert(nd->getAdjListIn().begin(), nd->getAdjListIn().end());
      }
    }
  }

  if (shrdNd == edg2->getFrom()) {
    if (_handlesA.count(edg2)) {
      for (auto nd : _handlesA.find(edg2)->second) {
        to.insert(nd->getAdjListIn().begin(), nd->getAdjListIn().end());
      }
    }
  } else {
    if (_handlesB.count(edg2)) {
      for (auto nd : _handlesB.find(edg2)->second) {
        to.insert(nd->getAdjListIn().begin(), nd->getAdjListIn().end());
      }
    }
  }

  EDijkstra::EList<RestrNodePL, RestrEdgePL> resEdges;

  CostFunc cFunc(r, 100000);
  double cost = EDijkstra::shortestPath(
      from, to, cFunc,
      EDijkstra::ZeroHeurFunc<RestrNodePL, RestrEdgePL, double>(), &resEdges);

  std::cerr << "From " << edg1 << " to " << edg2 << " and line "
            << r->getLabel() << "(" << r->getId() << "), " << from.size()
            << " x " << to.size() << " orig distance is ";
  std::cerr << cost << ", cur distance is " << curDist << std::endl;
  std::cerr << "From: ";
  for (auto edg : from) std::cerr << edg << ", ";
  std::cerr << std::endl << "To: ";
  for (auto edg : to) std::cerr << edg << ", ";
  std::cerr << std::endl;

  std::cerr << "Edges: " << std::endl;
  for (auto edg : resEdges) {
    std::cerr << edg << std::endl;
  }

  return cost / curDist < 1.5;
}
