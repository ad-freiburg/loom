// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <climits>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/statinserter/StatInserter.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace topo;

using topo::config::TopoConfig;

using util::geo::Point;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Box;
using util::geo::DBox;
using util::geo::extendBox;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::Station;
using shared::transitgraph::TransitEdge;
using shared::transitgraph::TransitEdgePair;
using shared::transitgraph::TransitNodePL;
using shared::transitgraph::TransitEdgePL;

// _____________________________________________________________________________
StatInserter::StatInserter(const TopoConfig* cfg, TransitGraph* g)
    : _cfg(cfg), _g(g) {}

// _____________________________________________________________________________
EdgeGrid StatInserter::geoIndex() {
  EdgeGrid grid(120, 120, bbox());

  for (auto n : *_g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      grid.add(e->pl().getPolyline().getLine(), e);
    }
  }

  return grid;
}

// _____________________________________________________________________________
DBox StatInserter::bbox() const {
  DBox b;

  for (auto n : *_g->getNds()) {
    b = extendBox(*n->pl().getGeom(), b);
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      b = extendBox(e->pl().getPolyline().getLine(), b);
    }
  }

  return b;
}

// _____________________________________________________________________________
void StatInserter::init() {
  for (auto nd : *_g->getNds()) {
    if (nd->pl().getStops().size()) {
      StationOcc occ{nd->pl().getStops().front(),
                     {nd->getAdjList().begin(), nd->getAdjList().end()}};
      _statClusters.push_back({occ});
      nd->pl().clearStops();
    }
  }
  std::cerr << "Collected " << _statClusters.size() << " station clusters..."
            << std::endl;
}

// _____________________________________________________________________________
std::pair<size_t, size_t> StatInserter::served(
    const std::vector<TransitEdge*>& adj,
    const std::set<const TransitEdge*>& toServe, const OrigEdgs& origEdgs) {
  std::set<const TransitEdge*> contained;

  for (auto e : adj)
    contained.insert(origEdgs.find(e)->second.begin(),
                     origEdgs.find(e)->second.end());

  std::set<const TransitEdge*> iSect;
  set_intersection(contained.begin(), contained.end(), toServe.begin(),
                   toServe.end(), std::inserter(iSect, iSect.begin()));

  // TODO: second return value should be falsely served lines
  return {iSect.size(), 0};
}

// _____________________________________________________________________________
double StatInserter::candScore(const StationCand& c) {
  double score = 0;

  score += c.dist;
  score += static_cast<double>(c.shouldServ) /
           static_cast<double>(c.truelyServ) * 100;

  // add a penalty if a station is too close to an existing node
  if (c.edg && c.edg->pl().getPolyline().getLength() * (1 - c.pos) < 100)
    score += 200;
  if (c.edg && c.edg->pl().getPolyline().getLength() * (c.pos) < 100)
    score += 200;

  return score;
}

// _____________________________________________________________________________
std::vector<StationCand> StatInserter::candidates(const StationOcc& occ,
                                                  const EdgeGrid& idx,
                                                  const OrigEdgs& origEdgs) {
  std::vector<StationCand> ret;
  std::set<TransitEdge*> neighbors;
  idx.get(util::geo::pad(util::geo::getBoundingBox(occ.station.pos), 100),
          &neighbors);

  std::cerr << "Got " << neighbors.size() << " candidates..." << std::endl;

  for (auto edg : neighbors) {
    auto pos = edg->pl().getPolyline().projectOn(occ.station.pos);
    double d = util::geo::dist(pos.p, occ.station.pos);

    size_t truelyServed, falselyServed;
    std::tie(truelyServed, falselyServed) = served({edg}, occ.edges, origEdgs);

    ret.push_back(StationCand{edg, pos.totalPos, 0, d, occ.edges.size(),
                              truelyServed, falselyServed});

    std::tie(truelyServed, falselyServed) =
        served(edg->getFrom()->getAdjList(), occ.edges, origEdgs);
    ret.push_back(StationCand{
        0, 0, edg->getFrom(),
        util::geo::dist(*edg->getFrom()->pl().getGeom(), occ.station.pos),
        occ.edges.size(), truelyServed, falselyServed});
    std::tie(truelyServed, falselyServed) =
        served(edg->getTo()->getAdjList(), occ.edges, origEdgs);
    ret.push_back(StationCand{
        0, 0, edg->getTo(),
        util::geo::dist(*edg->getTo()->pl().getGeom(), occ.station.pos),
        occ.edges.size(), truelyServed, falselyServed});
  }

  struct {
    bool operator()(const StationCand& a, const StationCand& b) {
      return candScore(a) < candScore(b);
    }
  } cmp;

  std::sort(ret.begin(), ret.end(), cmp);

  std::cerr << "  Cands: " << std::endl;
  for (auto cand : ret) {
    if (cand.edg) {
      std::cerr << "    Edg " << cand.edg << " at position " << cand.pos
                << " with dist = " << cand.dist << " truely serving "
                << cand.truelyServ << " edges, falsely serving "
                << cand.falselyServ << " edges." << std::endl;
    } else {
      std::cerr << "    Nd " << cand.nd << " with dist = " << cand.dist
                << " truely serving " << cand.truelyServ
                << " edges, falsely serving " << cand.falselyServ << " edges."
                << std::endl;
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool StatInserter::insertStations(const OrigEdgs& origEdgs) {
  auto idx = geoIndex();

  for (auto st : _statClusters) {
    if (st.size() == 0) continue;

    auto repr = st.front();
    std::cerr << "Inserting " << repr.station.name << std::endl;

    auto cands = candidates(repr, idx, origEdgs);

    if (cands.size() == 0) continue;

    if (cands.front().edg) {
      auto e = cands.front().edg;

      auto spl = split(e->pl(), e->getFrom(), e->getTo(), cands.front().pos);

      shared::transitgraph::TransitGraph::sharedNode(spl.first, spl.second)
          ->pl()
          .addStop(repr.station);

      idx.add(*spl.first->pl().getGeom(), spl.first);
      idx.add(*spl.second->pl().getGeom(), spl.second);

      edgeRpl(e->getFrom(), e, spl.first);
      edgeRpl(e->getTo(), e, spl.second);

      _g->delEdg(e->getFrom(), e->getTo());
      idx.remove(e);
    } else {
      cands.front().nd->pl().addStop(repr.station);
    }
  }

  return true;
}

// _____________________________________________________________________________
TransitEdgePair StatInserter::split(TransitEdgePL& a, TransitNode* fr,
                                    TransitNode* to, double p) {
  TransitEdge* ret;
  auto right = a.getPolyline().getSegment(p, 1);
  a.setPolyline(a.getPolyline().getSegment(0, p));
  auto helper = _g->addNd(a.getPolyline().back());
  auto ro = a.getRoutes().begin();
  auto helperEdg = _g->addEdg(helper, to, right);

  while (ro != a.getRoutes().end()) {
    if (ro->direction == to) {
      auto* route = ro->route;  // store because of deletion below
      ro = a.getRoutes().erase(ro);
      a.addRoute(route, helper);
      helperEdg->pl().addRoute(route, to);
    } else if (ro->direction == fr) {
      helperEdg->pl().addRoute(ro->route, helper);
      ro++;
    } else {
      helperEdg->pl().addRoute(ro->route, 0);
      ro++;
    }
  }

  ret = _g->addEdg(fr, helper, a);

  return {ret, helperEdg};
}

// _____________________________________________________________________________
void StatInserter::edgeRpl(TransitNode* n, const TransitEdge* oldE,
                           const TransitEdge* newE) {
  if (oldE == newE) return;
  // replace in from
  for (auto& r : n->pl().getConnExc()) {
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
  for (auto& r : n->pl().getConnExc()) {
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
