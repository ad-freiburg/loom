// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include <unordered_map>
#include "shared/linegraph/LineGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/restr/RestrGraph.h"
#include "topo/restr/RestrInferrer.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Dijkstra.h"

using topo::restr::RestrInferrer;

// _____________________________________________________________________________
RestrInferrer::RestrInferrer(const TopoConfig* cfg, LineGraph* g)
    : _cfg(cfg), _tg(g) {}

// _____________________________________________________________________________
void RestrInferrer::init() {
  for (auto nd : *_tg->getNds()) {
    _nMap[nd] = _rg.addNd();
  }

  for (auto nd : *_tg->getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      const auto& pl = edg->pl().getPolyline();
      _eMap[edg] = {_rg.addEdg(_nMap[edg->getFrom()], _nMap[edg->getTo()], pl),
                    _rg.addEdg(_nMap[edg->getTo()], _nMap[edg->getFrom()],
                               pl.reversed())};

      for (auto r : edg->pl().getRoutes()) {
        if (r.direction == 0 || r.direction == edg->getTo()) {
          _eMap[edg][0]->pl().routes.insert(r.route);
        }

        if (r.direction == 0 || r.direction == edg->getFrom()) {
          _eMap[edg][1]->pl().routes.insert(r.route);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void RestrInferrer::infer(const OrigEdgs& origEdgs) {
  addHndls(origEdgs);

  for (auto nd : *_tg->getNds()) {
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

          if (!check(ro1.route, edg1, edg2) && !check(ro1.route, edg2, edg1)) {
            nd->pl().addConnExc(ro1.route, edg1, edg2);
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void RestrInferrer::addHndls(const OrigEdgs& origEdgs) {
  std::map<RestrEdge*, HndlLst> handles;

  // collect the handles
  for (auto nd : *_tg->getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      addHndls(edg, origEdgs, &handles);
    }
  }

  // add them to the edge
  for (auto edgHndl : handles) {
    std::sort(edgHndl.second.begin(), edgHndl.second.end(), HndlCmp());

    auto lastNd = edgHndl.first->getFrom();
    double lastPos = 0;
    for (auto hndl : edgHndl.second) {
      auto e =
          _rg.addEdg(lastNd, hndl.first,
                     edgHndl.first->pl().geom.getSegment(lastPos, hndl.second));
      e->pl().routes = edgHndl.first->pl().routes;

      lastNd = hndl.first;
      lastPos = hndl.second;
    }
    auto e = _rg.addEdg(lastNd, edgHndl.first->getTo(),
                        edgHndl.first->pl().geom.getSegment(lastPos, 1));
    e->pl().routes = edgHndl.first->pl().routes;

    _rg.delEdg(edgHndl.first->getFrom(), edgHndl.first->getTo());
  }
}

// _____________________________________________________________________________
void RestrInferrer::addHndls(const LineEdge* e, const OrigEdgs& origEdgs,
                             std::map<RestrEdge*, HndlLst>* handles) {
  double MAX_DIST = _cfg->maxAggrDistance;

  auto hndlPA = e->pl().getPolyline().getPointAt(1.0 / 3.0).p;
  auto hndlPB = e->pl().getPolyline().getPointAt(2.0 / 3.0).p;

  auto hndlLA =
      e->pl()
          .getPolyline()
          .getOrthoLineAtDist(MAX_DIST, e->pl().getRoutes().size() * MAX_DIST)
          .getLine();
  auto hndlLB =
      e->pl()
          .getPolyline()
          .getOrthoLineAtDist(e->pl().getPolyline().getLength() - MAX_DIST,
                              e->pl().getRoutes().size() * MAX_DIST)
          .getLine();

  for (auto edg : origEdgs.find(e)->second) {
    auto origFr = const_cast<LineEdge*>(edg);

    const auto& edgs = _eMap.find(origFr)->second;

    for (auto restrE : edgs) {
      if (util::geo::intersects(hndlLA, restrE->pl().geom.getLine())) {
        auto projA = restrE->pl().geom.projectOn(hndlPA);
        auto handleNdA = _rg.addNd();

        _handlesA[e].insert(handleNdA);
        (*handles)[restrE].push_back({handleNdA, projA.totalPos});
      }

      if (util::geo::intersects(hndlLB, restrE->pl().geom.getLine())) {
        auto projB = restrE->pl().geom.projectOn(hndlPB);
        auto handleNdB = _rg.addNd();

        _handlesB[e].insert(handleNdB);
        (*handles)[restrE].push_back({handleNdB, projB.totalPos});
      }
    }
  }
}

// _____________________________________________________________________________
bool RestrInferrer::check(const Route* r, const LineEdge* edg1,
                          const LineEdge* edg2) const {
  std::set<RestrEdge *> from, to;
  auto shrdNd = shared::linegraph::LineGraph::sharedNode(edg1, edg2);

  double curD = edg1->pl().getPolyline().getLength() * 0.33 +
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

  // curdist + 500 is the inf. We do not have to check any further as we
  // only return true below if cost - curD < 500 <=> cost < curD + 500
  CostFunc cFunc(r, curD + _cfg->maxLengthDev);
  return EDijkstra::shortestPath(from, to, cFunc) - curD < _cfg->maxLengthDev;
}
