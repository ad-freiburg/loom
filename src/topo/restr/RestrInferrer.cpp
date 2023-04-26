// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include <unordered_map>
#include "shared/linegraph/LineGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/mapconstructor/MapConstructor.h"
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
  for (auto nd : _tg->getNds()) {
    _nMap[nd] = _rg.addNd(*nd->pl().getGeom());
  }

  for (auto nd : _tg->getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      const auto& pl = edg->pl().getPolyline();
      _eMap[edg] = {_rg.addEdg(_nMap[edg->getFrom()], _nMap[edg->getTo()], pl),
                    _rg.addEdg(_nMap[edg->getTo()], _nMap[edg->getFrom()],
                               pl.reversed())};

      for (auto r : edg->pl().getLines()) {
        if (r.direction == 0 || r.direction == edg->getTo()) {
          _eMap[edg][0]->pl().lines.insert(r.line);
        }

        if (r.direction == 0 || r.direction == edg->getFrom()) {
          _eMap[edg][1]->pl().lines.insert(r.line);
        }
      }
    }
  }

  // copy turn restrictions from original graph
  for (auto nd : _tg->getNds()) {
    for (auto ex : nd->pl().getConnExc()) {
      auto line = ex.first;
      for (auto exPair : ex.second) {
        const LineEdge* edgeFr = exPair.first;

        for (RestrEdge* rEdgeFr : _eMap[edgeFr]) {
          for (auto edgeTo : exPair.second) {
            for (RestrEdge* rEdgeTo : _eMap[edgeTo]) {
              _nMap[nd]->pl().restrs[line][rEdgeFr].insert(rEdgeTo);
              _nMap[nd]->pl().restrs[line][rEdgeTo].insert(rEdgeFr);
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
size_t RestrInferrer::infer(const OrigEdgs& origEdgs) {
  // delete all existing restrictions

  for (auto nd : _tg->getNds()) {
    nd->pl().clearConnExc();
  }

  addHndls(origEdgs);

  size_t ret = 0;

  // debug output
  // util::geo::output::GeoGraphJsonOutput out;
  // std::ofstream outs;
  // outs.open("restr_graph.json");
  // out.printLatLng(_rg, outs);

  for (auto nd : _tg->getNds()) {
    for (auto edg1 : nd->getAdjList()) {
      // check every other edge
      for (auto edg2 : nd->getAdjList()) {
        if (edg1 == edg2) continue;

        for (auto ro1 : edg1->pl().getLines()) {
          if (!edg2->pl().hasLine(ro1.line)) continue;

          const auto& ro2 = edg2->pl().lineOcc(ro1.line);

          if (ro1.direction != 0 && ro2.direction != 0 &&
              ro1.direction == ro2.direction)
            continue;

          if (ro1.direction != 0 && ro2.direction != 0 &&
              edg1->getOtherNd(ro1.direction) ==
                  edg2->getOtherNd(ro2.direction)) {
            continue;
          }

          if (!check(ro1.line, edg1, edg2) && !check(ro1.line, edg2, edg1)) {
            nd->pl().addConnExc(ro1.line, edg1, edg2);
            ret++;
          }
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
void RestrInferrer::addHndls(const OrigEdgs& origEdgs) {
  std::map<RestrEdge*, HndlLst> handles;

  // collect the handles
  for (auto nd : _tg->getNds()) {
    for (auto edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      addHndls(edg, origEdgs, &handles);
    }
  }

  // add them to the edge
  for (auto edgHndl : handles) {
    // sort the handles by their occurance on the line edge
    std::sort(edgHndl.second.begin(), edgHndl.second.end(), HndlCmp());

    auto lastNd = edgHndl.first->getFrom();
    double lastPos = 0;

    for (auto hndl : edgHndl.second) {
      auto e =
          _rg.addEdg(lastNd, hndl.first,
                     edgHndl.first->pl().geom.getSegment(lastPos, hndl.second));
      e->pl().lines = edgHndl.first->pl().lines;

      lastNd = hndl.first;
      lastPos = hndl.second;

      // replace any existing exception occurances of the old edge with the new
      edgeRpl(e->getFrom(), edgHndl.first, e);
      edgeRpl(e->getTo(), edgHndl.first, e);
    }

    // last edge
    auto e = _rg.addEdg(lastNd, edgHndl.first->getTo(),
                        edgHndl.first->pl().geom.getSegment(lastPos, 1));
    e->pl().lines = edgHndl.first->pl().lines;

    // replace any existing exception occurances of the old edge with the new
    edgeRpl(e->getFrom(), edgHndl.first, e);
    edgeRpl(e->getTo(), edgHndl.first, e);

    // delete original edge
    _rg.delEdg(edgHndl.first->getFrom(), edgHndl.first->getTo());
  }
}

// _____________________________________________________________________________
void RestrInferrer::edgeRpl(RestrNode* n, const RestrEdge* oldE,
                            const RestrEdge* newE) {
  if (oldE == newE) return;
  // replace in from
  for (auto& r : n->pl().restrs) {
    auto exFr = r.second.begin();
    while (exFr != r.second.end()) {
      if (exFr->first == oldE) {
        std::swap(r.second[newE], exFr->second);
        exFr = r.second.erase(exFr);
      } else {
        exFr++;
      }
    }
  }

  // replace in to
  for (auto& r : n->pl().restrs) {
    for (auto& exFr : r.second) {
      auto exTo = exFr.second.begin();
      while (exTo != exFr.second.end()) {
        if (*exTo == oldE) {
          exFr.second.insert(newE);
          exTo = exFr.second.erase(exTo);
        } else {
          exTo++;
        }
      }
    }
  }
}

// _____________________________________________________________________________
void RestrInferrer::addHndls(const LineEdge* e, const OrigEdgs& origEdgs,
                             std::map<RestrEdge*, HndlLst>* handles) {
  // double MAX_DIST = _cfg->maxAggrDistance;
  AggrDistFunc aggrD(_cfg->maxAggrDistance);

  // auto hndlPA = e->pl().getPolyline().getPointAt(1.0 / 3.0).p;
  // auto hndlPB = e->pl().getPolyline().getPointAt(2.0 / 3.0).p;

  // we add a small buffer to account for the skip heuristic in the
  // shared segments collapsing
  double MAX_DIST = _cfg->maxTurnRestrCheckDist * 2;
  double checkPos =
      std::min(e->pl().getPolyline().getLength() / 2, 2 * _cfg->maxAggrDistance);

  auto hndlLA =
      e->pl().getPolyline().getOrthoLineAt(1.0 / 3.0, MAX_DIST).getLine();
  auto hndlLB =
      e->pl().getPolyline().getOrthoLineAt(2.0 / 3.0, MAX_DIST).getLine();

  auto hndlLACheck =
      e->pl().getPolyline().getOrthoLineAtDist(checkPos, MAX_DIST).getLine();
  auto hndlLBCheck = e->pl()
                    .getPolyline()
                    .getOrthoLineAtDist(
                        e->pl().getPolyline().getLength() - checkPos, MAX_DIST)
                    .getLine();

  auto a =_rg.addNd(hndlLA.front());
  auto b = _rg.addNd(hndlLA.back());
  _rg.addEdg(a, b, RestrEdgePL(hndlLA));

  a =_rg.addNd(hndlLB.front());
  b = _rg.addNd(hndlLB.back());
  _rg.addEdg(a, b, RestrEdgePL(hndlLB));

   // a =_rg.addNd(hndlLACheck.front());
   // b = _rg.addNd(hndlLACheck.back());
  // _rg.addEdg(a, b, RestrEdgePL(hndlLACheck));

  // a =_rg.addNd(hndlLBCheck.front());
  // b = _rg.addNd(hndlLBCheck.back());
  // _rg.addEdg(a, b, RestrEdgePL(hndlLBCheck));

  for (auto edg : origEdgs.find(e)->second) {
    auto origFr = const_cast<LineEdge*>(edg);
    const auto& edgs = _eMap.find(origFr)->second;

    assert(edgs.size());

    for (auto restrE : edgs) {
      // the geometry of the original line
      auto geom = restrE->pl().geom.getLine();

      // make sure this geometry is connected to the nodes
      geom.insert(geom.begin(), *restrE->getFrom()->pl().getGeom());
      geom.insert(geom.end(), *restrE->getTo()->pl().getGeom());

      if (util::geo::intersects(hndlLACheck, geom)) {
        auto isects = util::geo::intersection(geom, hndlLA);

        // if no intersections found, fall back to the handlLACheck
        if (isects.size() == 0) isects = util::geo::intersection(geom, hndlLACheck);

        for (const auto& isect : isects) {
          auto projA = restrE->pl().geom.projectOn(isect);
          auto handleNdA = _rg.addNd(projA.p);

          _handlesA[e].insert(handleNdA);
          (*handles)[restrE].push_back({handleNdA, projA.totalPos});
        }
      }

      if (util::geo::intersects(hndlLBCheck, geom)) {
        auto isects = util::geo::intersection(geom, hndlLB);

        // if no intersections found, fall back to the handlLACheck
        if (isects.size() == 0) isects = util::geo::intersection(geom, hndlLBCheck);

        for (const auto& isect : isects) {
          auto projB = restrE->pl().geom.projectOn(isect);
          auto handleNdB = _rg.addNd(projB.p);

          _handlesB[e].insert(handleNdB);
          (*handles)[restrE].push_back({handleNdB, projB.totalPos});
        }
      }
    }
  }
}

// _____________________________________________________________________________
bool RestrInferrer::check(const Line* r, const LineEdge* edg1,
                          const LineEdge* edg2) const {
  std::set<RestrEdge*> from, to;
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

  // curdist + maxL is the inf. We do not have to check any further as we
  // only return true below if cost - curD < maxL <=> cost < curD + maxL
  // + epsilon to avoid integer rounding issues in the < comparison below
  double eps = 0.1;
  CostFunc cFunc(r, curD + _cfg->maxLengthDev + eps, _cfg->turnInferFullTurnPen);

  return EDijkstra::shortestPath(from, to, cFunc) - curD < _cfg->maxLengthDev;
}
