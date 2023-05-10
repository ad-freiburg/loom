// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <climits>

#include "shared/linegraph/LineGraph.h"
#include "topo/statinserter/StatInserter.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace topo;

using topo::config::TopoConfig;

using util::geo::Box;
using util::geo::DBox;
using util::geo::DPoint;
using util::geo::extendBox;
using util::geo::Grid;
using util::geo::Point;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePair;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::Station;

// _____________________________________________________________________________
StatInserter::StatInserter(const TopoConfig* cfg, LineGraph* g)
    : _cfg(cfg), _g(g) {
  UNUSED(_cfg);
}

// _____________________________________________________________________________
EdgeGeoIdx StatInserter::geoIndex() {
  EdgeGeoIdx grid;

  for (auto n : _g->getNds()) {
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

  for (auto n : _g->getNds()) {
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
  std::map<std::string, size_t> existing;

  for (auto nd : _g->getNds()) {
    if (nd->pl().stops().size()) {
      auto stop = nd->pl().stops().front();

      auto exI = existing.find(stop.name);
      if (exI != existing.end()) {
        auto& ex = _statClusters[exI->second];
        ex.front().edges.insert(nd->getAdjList().begin(),
                                nd->getAdjList().end());
        for (auto e : nd->getAdjList()) {
          for (const auto& lo : e->pl().getLines()) {
            ex.front().lines.insert(lo.line);
          }
        }
      } else {
        StationOcc occ{
            stop, {nd->getAdjList().begin(), nd->getAdjList().end()}, {}};
        for (auto e : nd->getAdjList()) {
          for (const auto& lo : e->pl().getLines()) {
            occ.lines.insert(lo.line);
          }
        }
        _statClusters.push_back({occ});
        existing[stop.name] = _statClusters.size() - 1;
      }

      nd->pl().clearStops();
    }
  }
  LOGTO(DEBUG, std::cerr) << "Collected " << _statClusters.size()
                          << " station clusters...";
}

// _____________________________________________________________________________
StationOcc StatInserter::unserved(const std::vector<LineEdge*>& adj,
                                  const StationOcc& stationOcc,
                                  const OrigEdgs& origEdgs) {
  StationOcc ret{stationOcc.station, {}, {}};
  std::set<const LineEdge*> contained;
  std::set<const shared::linegraph::Line*> containedLines;

  for (auto e : adj)
    contained.insert(origEdgs.find(e)->second.begin(),
                     origEdgs.find(e)->second.end());

  for (auto e : adj) {
    for (auto lo : e->pl().getLines()) {
      containedLines.insert(lo.line);
    }
  }

  std::set<const LineEdge*> diff;
  set_difference(stationOcc.edges.begin(), stationOcc.edges.end(),
                 contained.begin(), contained.end(),
                 std::inserter(diff, diff.begin()));

  std::set<const shared::linegraph::Line*> diffLines;
  set_difference(stationOcc.lines.begin(), stationOcc.lines.end(),
                 containedLines.begin(), containedLines.end(),
                 std::inserter(diffLines, diffLines.begin()));

  ret.edges.insert(diff.begin(), diff.end());
  ret.lines.insert(diffLines.begin(), diffLines.end());

  return ret;
}

// _____________________________________________________________________________
std::pair<size_t, size_t> StatInserter::served(
    const std::vector<LineEdge*>& adj, const std::set<const LineEdge*>& toServe,
    const std::set<const shared::linegraph::Line*>& linesToServe,
    const OrigEdgs& origEdgs) {
  std::set<const LineEdge*> contained;
  std::set<const shared::linegraph::Line*> containedLines;

  for (auto e : adj)
    contained.insert(origEdgs.find(e)->second.begin(),
                     origEdgs.find(e)->second.end());

  for (auto e : adj) {
    for (auto lo : e->pl().getLines()) {
      containedLines.insert(lo.line);
    }
  }

  std::set<const LineEdge*> iSect;
  set_intersection(contained.begin(), contained.end(), toServe.begin(),
                   toServe.end(), std::inserter(iSect, iSect.begin()));

  std::set<const shared::linegraph::Line*> iSectLines;
  set_intersection(containedLines.begin(), containedLines.end(),
                   linesToServe.begin(), linesToServe.end(),
                   std::inserter(iSectLines, iSectLines.begin()));

  return {iSect.size(), iSectLines.size()};
}

// _____________________________________________________________________________
double StatInserter::candScore(const StationCand& c) {
  double score = 0;

  score += c.dist;
  score += static_cast<double>(c.shouldServ) /
           static_cast<double>(c.truelyServ) * 100;
  score += static_cast<double>(c.shouldServLines) /
           static_cast<double>(c.truelyServedLines) * 200;

  // add a penalty if a station is inserted to an edge, and is
  // too close to an existing node (prefer the existing node later on)
  if (c.edg && c.edg->pl().getPolyline().getLength() * (1 - c.pos) < 100)
    score += 200;
  if (c.edg && c.edg->pl().getPolyline().getLength() * (c.pos) < 100)
    score += 200;

  return score;
}

// _____________________________________________________________________________
std::vector<StationCand> StatInserter::candidates(const StationOcc& occ,
                                                  const EdgeGeoIdx& idx,
                                                  const OrigEdgs& origEdgs) {
  std::vector<StationCand> ret;
  std::set<LineEdge*> neighbors;
  idx.get(util::geo::pad(util::geo::getBoundingBox(occ.station.pos), 250),
          &neighbors);

  LOGTO(VDEBUG, std::cerr) << "Got " << neighbors.size() << " candidates...";

  for (auto edg : neighbors) {
    auto pos = edg->pl().getPolyline().projectOn(occ.station.pos);
    double d = util::geo::dist(pos.p, occ.station.pos);

    size_t truelyServed, truelyServedLines;

    // add whole edge as cand
    std::tie(truelyServed, truelyServedLines) =
        served({edg}, occ.edges, occ.lines, origEdgs);
    StationOcc remaining = unserved({edg}, occ, origEdgs);
    ret.push_back(StationCand{edg, pos.totalPos, 0, d, occ.edges.size(),
                              occ.lines.size(), truelyServed, truelyServedLines,
                              remaining});

    // add from node as cand
    std::tie(truelyServed, truelyServedLines) =
        served(edg->getFrom()->getAdjList(), occ.edges, occ.lines, origEdgs);
    remaining = unserved(edg->getFrom()->getAdjList(), occ, origEdgs);
    ret.push_back(StationCand{
        0, 0, edg->getFrom(),
        util::geo::dist(*edg->getFrom()->pl().getGeom(), occ.station.pos),
        occ.edges.size(), occ.lines.size(), truelyServed, truelyServedLines,
        remaining});

    // add to node as cand
    std::tie(truelyServed, truelyServedLines) =
        served(edg->getTo()->getAdjList(), occ.edges, occ.lines, origEdgs);
    remaining = unserved(edg->getTo()->getAdjList(), occ, origEdgs);
    ret.push_back(StationCand{
        0, 0, edg->getTo(),
        util::geo::dist(*edg->getTo()->pl().getGeom(), occ.station.pos),
        occ.edges.size(), occ.lines.size(), truelyServed, truelyServedLines,
        remaining});
  }

  struct {
    bool operator()(const StationCand& a, const StationCand& b) {
      return candScore(a) < candScore(b);
    }
  } cmp;

  std::sort(ret.begin(), ret.end(), cmp);

  LOGTO(VDEBUG, std::cerr) << "  Cands: ";
  for (auto cand : ret) {
    if (cand.edg) {
      LOGTO(VDEBUG, std::cerr)
          << "    Edg " << cand.edg << " at position " << cand.pos
          << " with dist = " << cand.dist << " truely serving "
          << cand.truelyServ << "/" << occ.edges.size()
          << " edges, truely serving " << cand.truelyServedLines << "/"
          << occ.lines.size() << " lines.";
    } else {
      LOGTO(VDEBUG, std::cerr)
          << "    Nd " << cand.nd << " with dist = " << cand.dist
          << " truely serving " << cand.truelyServ << "/" << occ.edges.size()
          << " edges, truely serving " << cand.truelyServedLines << "/"
          << occ.lines.size() << " lines.";
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool StatInserter::insertStations(const OrigEdgs& origEdgs) {
  OrigEdgs modOrigEdgs = origEdgs;
  auto idx = geoIndex();

  for (auto st : _statClusters) {
    if (st.size() == 0) continue;

    auto curOcc = st.front();
    LOGTO(DEBUG, std::cerr) << "Inserting " << curOcc.station.name;

    int MAX_INSERTS = 3;
    int i = 0;

    while (i++ < MAX_INSERTS) {
      auto cands = candidates(curOcc, idx, modOrigEdgs);

      if (cands.size() == 0) {
        LOGTO(DEBUG, std::cerr) << "  (No insertion candidate found.)";
        break;
      }

      auto curCan = cands.front();

      if (curCan.truelyServ == 0) {
        LOGTO(DEBUG, std::cerr) << "  (No insertion candidate found.)";
        break;
      }

      if (curCan.edg) {
        auto e = curCan.edg;

        auto spl = split(e->pl(), e->getFrom(), e->getTo(), curCan.pos);

        shared::linegraph::LineGraph::sharedNode(spl.first, spl.second)
            ->pl()
            .addStop(curOcc.station);

        idx.add(*spl.first->pl().getGeom(), spl.first);
        idx.add(*spl.second->pl().getGeom(), spl.second);

        // UPDATE ORIGEDGES
        modOrigEdgs[spl.first] = modOrigEdgs[e];
        modOrigEdgs[spl.second] = modOrigEdgs[e];

        edgeRpl(e->getFrom(), e, spl.first);
        edgeRpl(e->getTo(), e, spl.second);

        _g->delEdg(e->getFrom(), e->getTo());
        idx.remove(e);
      } else {
        curCan.nd->pl().addStop(curOcc.station);
      }

      if (curCan.unserved.edges.size() == 0 &&
          curCan.unserved.lines.size() == 0)
        break;

      LOGTO(DEBUG, std::cerr)
          << "  inserting for remaining " << curCan.unserved.edges.size()
          << " unserved edges and/or " << curCan.unserved.lines.size()
          << " unserved lines...";
      curOcc = curCan.unserved;
    }
  }

  return true;
}

// _____________________________________________________________________________
LineEdgePair StatInserter::split(LineEdgePL& a, LineNode* fr, LineNode* to,
                                 double p) {
  LineEdge* ret;
  auto right = a.getPolyline().getSegment(p, 1);
  a.setPolyline(a.getPolyline().getSegment(0, p));
  auto helper = _g->addNd(a.getPolyline().back());
  auto helperEdg = _g->addEdg(helper, to, right);

  for (size_t i = 0; i < a.getLines().size(); i++) {
    auto ro = a.getLines()[i];
    if (ro.direction == to) {
      auto* route = ro.line;  // store because of deletion below
      a.delLine(ro.line);
      a.addLine(route, helper);
      helperEdg->pl().addLine(route, to);
      i--;
    } else if (ro.direction == fr) {
      helperEdg->pl().addLine(ro.line, helper);
    } else {
      helperEdg->pl().addLine(ro.line, 0);
    }
  }

  ret = _g->addEdg(fr, helper, a);

  return {ret, helperEdg};
}

// _____________________________________________________________________________
void StatInserter::edgeRpl(LineNode* n, const LineEdge* oldE,
                           const LineEdge* newE) {
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
