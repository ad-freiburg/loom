// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <unordered_map>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/restr/RestrGraph.h"
#include "topo/restr/RestrInferrer.h"
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
    }
  }
}

// _____________________________________________________________________________
void RestrInferrer::infer(TransitGraph* g, const OrigEdgs& origEdgs) {
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
                  edg2->getOtherNd(ro2.direction))
            continue;

          if (!check(ro1.route, edg1, edg2, origEdgs)) {
            nd->pl().addConnExc(ro1.route, edg1, edg2);
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
bool RestrInferrer::check(const Route* r, const TransitEdge* edg1,
                          const TransitEdge* edg2,
                          const OrigEdgs& origEdgs) const {
  // TODO: use max_aggr_dist here
  double MAX_DIST = 30;

  std::set<RestrEdge*> from;
  std::set<RestrEdge*> to;

  std::map<const RestrEdge*, double> sourcePos;
  std::map<const RestrEdge*, double> targetPos;

  auto shrdNd = shared::transitgraph::TransitGraph::sharedNode(edg1, edg2);

  util::geo::DPoint posFr, posTo;

  if (shrdNd == edg1->getFrom())
    posFr = edg1->pl().getPolyline().getPointAt(0.33).p;
  else
    posFr = edg1->pl().getPolyline().getPointAt(0.66).p;

  if (shrdNd == edg2->getFrom())
    posTo = edg2->pl().getPolyline().getPointAt(0.33).p;
  else
    posTo = edg2->pl().getPolyline().getPointAt(0.66).p;

  for (auto e : origEdgs.find(edg1)->second) {
    auto origFr = const_cast<TransitEdge*>(e);
    assert(_eMap.count(origFr));
    // if (util::geo::dist(posFr, *_eMap.find(origFr)->second->pl().getGeom()) <
        // MAX_DIST * edg1->pl().getRoutes().size()) {
      // from.insert(_eMap.find(origFr)->second);
      // sourcePos[origFr] = _eMap.find(origFr)
                              // ->second->pl()
                              // .getPolyline()
                              // .projectOn(posFr)
                              // .totalPos;
    // }
  // }

  // for (auto e : origEdgs.find(edg2)->second) {
    // auto origTo = const_cast<TransitEdge*>(e);
    // assert(_eMap.count(origTo));
    // if (util::geo::dist(posTo, *_eMap.find(origTo)->second->pl().getGeom()) <
        // MAX_DIST * edg2->pl().getRoutes().size()) {
      // to.insert(_eMap.find(origTo)->second);
      // targetPos[origTo] = _eMap.find(origTo)
                              // ->second->pl()
                              // .getPolyline()
                              // .projectOn(posTo)
                              // .totalPos;
    // }
  }

  CostFunc cFunc(sourcePos, targetPos);
  double cost = EDijkstra::shortestPath(
      from, to, cFunc,
      EDijkstra::ZeroHeurFunc<RestrNodePL, RestrEdgePL, double>());

  std::cerr << cost << std::endl;

  return false;
}
