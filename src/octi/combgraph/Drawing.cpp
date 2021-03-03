#include <iostream>
#include "octi/basegraph/BaseGraph.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
#include "util/geo/BezierCurve.h"
#include "util/graph/Dijkstra.h"
#include "util/log/Log.h"

using octi::basegraph::BaseGraph;
using octi::basegraph::GridEdge;
using octi::basegraph::GridEdgePL;
using octi::basegraph::GridNode;
using octi::basegraph::GridNodePL;
using octi::combgraph::CombEdge;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::Drawing;
using octi::combgraph::Score;
using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using util::geo::BezierCurve;
using util::graph::Dijkstra;

// _____________________________________________________________________________
double Drawing::score() const { return _c; }

// _____________________________________________________________________________
void Drawing::setBaseGraph(const BaseGraph* gg) { _gg = gg; }

// _____________________________________________________________________________
Score Drawing::fullScore() const {
  Score ret{0, 0, 0, 0, 0};

  for (auto c : _ndReachCosts) ret.move += c.second;
  for (auto c : _ndBndCosts) ret.bend += c.second;
  for (auto c : _edgCosts) ret.hop += c.second;
  for (auto c : _springCosts) ret.dense += c.second;
  ret.full = _c;

  return ret;
}

// _____________________________________________________________________________
void Drawing::draw(CombEdge* ce, const GrEdgList& ges, bool rev) {
  if (_c == std::numeric_limits<double>::infinity()) _c = 0;
  if (_edgs.count(ce)) _edgs[ce].clear();

  int l = 0;

  for (size_t i = 0; i < ges.size(); i++) {
    auto ge = ges[i];

    if (rev) {
      _nds[ce->getFrom()] =
          ges.front()->getTo()->pl().getParent()->pl().getId();
      _nds[ce->getTo()] = ges.back()->getFrom()->pl().getParent()->pl().getId();
    } else {
      _nds[ce->getTo()] = ges.front()->getTo()->pl().getParent()->pl().getId();
      _nds[ce->getFrom()] =
          ges.back()->getFrom()->pl().getParent()->pl().getId();
    }

    // there are three kinds of cost contained in a result:
    //  a) node reach costs, which model the cost it takes to move a node
    //     away from its original position. They are only added to the
    //     first edge reaching this node, to prevent a weight for nodes with
    //     higher degree. These costs cannot be removed if a settled edge
    //     is removed, and they are exactly the cost of the sink edge
    //  b) node bend costs at input nodes, which are induced by a single edge
    //     and which can be removed if a settled edge is removed.
    //  c) edge bend costs at non-input nodes, these can be safely removed if
    //     a settled edge is unsettled
    //  d) edge costs, which can also be safely removed if a settled edge is
    //     unsettled

    _c += ge->pl().cost();

    if (i == 0) {
      if (!_ndReachCosts.count(rev ? ce->getFrom() : ce->getTo())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[rev ? ce->getFrom() : ce->getTo()] = ge->pl().cost();
        _ndBndCosts[rev ? ce->getFrom() : ce->getTo()] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[rev ? ce->getFrom() : ce->getTo()] += ge->pl().cost();
        // assert(!ges[i + 1]->pl().isSecondary());
      }
    } else if (i == ges.size() - 1) {
      if (!_ndReachCosts.count(rev ? ce->getTo() : ce->getFrom())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[rev ? ce->getTo() : ce->getFrom()] = ge->pl().cost();
        _ndBndCosts[rev ? ce->getTo() : ce->getFrom()] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[rev ? ce->getTo() : ce->getFrom()] += ge->pl().cost();
        // assert(!ges[i - 1]->pl().isSecondary());
      }
    } else {
      if (!ge->pl().isSecondary()) l++;
      _edgCosts[ce] += ge->pl().cost();
    }

    if (rev) {
      auto e = _gg->getEdg(ges[ges.size() - 1 - i]->getTo(),
                           ges[ges.size() - 1 - i]->getFrom());

      if (!e->pl().isSecondary())
        _edgs[ce].push_back(
            {e->getFrom()->pl().getId(), e->getTo()->pl().getId()});
    } else {
      if (!ges[i]->pl().isSecondary())
        _edgs[ce].push_back(
            {ges[i]->getFrom()->pl().getId(), ges[i]->getTo()->pl().getId()});
    }
  }

  // variables named as in the EuroVis paper
  int k = ce->pl().getChilds().size() - 1;

  double c = _gg->getPens().densityPen / (k);

  double F = c * (k + 1 - l);
  double E = 0.5 * c * (k + 1 - l) * (k + 1 - l);

  double pen = 0;
  if (F > 0) pen = E;

  _springCosts[ce] = pen;
  _c += _springCosts[ce];
}

// _____________________________________________________________________________
const GridNode* Drawing::getGrNd(const CombNode* cn) {
  return _gg->getGrNdById(_nds[cn]);
}

// _____________________________________________________________________________
void Drawing::getLineGraph(LineGraph* target) const {
  std::map<LineNode*, LineNode*> m;
  std::map<GridNode*, LineNode*> mm;

  std::vector<Segment> pathSegments;
  std::map<const CombEdge*, std::vector<std::pair<size_t, bool>>> combEdgSegs;

  std::map<const GridEdge*, size_t> grEdgToPathSeg;

  // settle grid nodes
  for (auto ndpair : _nds) {
    auto n = ndpair.first;
    for (auto f : n->getAdjListOut()) {
      if (f->getFrom() != n) continue;
      if (_edgs.find(f) == _edgs.end()) {
        LOGTO(WARN, std::cerr) << "Edge " << f << " was not drawn, skipping...";
        continue;
      }

      auto path = _edgs.find(f)->second;
      GridNode* from =
          _gg->getGrEdgById(path.back())->getFrom()->pl().getParent();
      GridNode* to = _gg->getGrEdgById(path.front())->getTo()->pl().getParent();

      if (mm.find(from) == mm.end()) {
        auto payload = n->pl().getParent()->pl();
        payload.setGeom(*from->pl().getGeom());
        mm[from] = target->addNd(payload);
        m[n->pl().getParent()] = mm[from];
      }

      if (mm.find(to) == mm.end()) {
        auto payload = f->getTo()->pl().getParent()->pl();
        payload.setGeom(*to->pl().getGeom());
        mm[to] = target->addNd(payload);
        m[f->getTo()->pl().getParent()] = mm[to];
      }
    }
  }

  // build segments per path
  for (auto ndpair : _nds) {
    auto n = ndpair.first;
    for (auto f : n->getAdjListOut()) {
      if (f->getFrom() != n) continue;
      if (_edgs.find(f) == _edgs.end()) {
        LOGTO(WARN, std::cerr) << "Edge " << f << " was not drawn, skipping...";
        continue;
      }

      auto path = _edgs.find(f)->second;

      std::set<CombEdge*> curResEdgs;

      for (size_t i = 0; i < path.size(); i++) {
        auto pathEdge = path[i];
        const GridEdge* grEdg = _gg->getGrEdgById(pathEdge);

        const GridEdge* grCounterEdg =
            _gg->getGrEdgById({pathEdge.second, pathEdge.first});

        if (grEdgToPathSeg.count(grEdg)) {
          size_t existSegment = grEdgToPathSeg[grEdg];
          if (combEdgSegs[f].size() == 0 ||
              combEdgSegs[f].back().first != existSegment) {
            combEdgSegs[f].push_back({existSegment, false});
            continue;
          }
        } else if (grEdgToPathSeg.count(grCounterEdg)) {
          size_t existSegment = grEdgToPathSeg[grCounterEdg];
          if (combEdgSegs[f].size() == 0 ||
              combEdgSegs[f].back().first != existSegment) {
            combEdgSegs[f].push_back({existSegment, true});
            continue;
          }
        }

        auto newResEdgs = _gg->getResEdgsDirInd(grEdg);

        if (newResEdgs == curResEdgs) {
          pathSegments[combEdgSegs[f].back().first].start =
              grEdg->getFrom()->pl().getParent();
          pathSegments[combEdgSegs[f].back().first].path.push_back(pathEdge);
          assert(pathSegments.back().end != pathSegments.back().start);
        } else {
          pathSegments.push_back({grEdg->getFrom()->pl().getParent(),
                                  grEdg->getTo()->pl().getParent(),
                                  {pathEdge},
                                  {},
                                  {},
                                  newResEdgs});
          assert(pathSegments.back().end != pathSegments.back().start);
          combEdgSegs[f].push_back({pathSegments.size() - 1, false});
          curResEdgs = newResEdgs;
        }
        grEdgToPathSeg[grEdg] = combEdgSegs[f].back().first;
      }
    }
  }

  // write geometries of segments
  for (auto& segment : pathSegments)
    segment.geom = _gg->geomFromPath(segment.path);

  // add nodes to segments
  for (auto ndpair : _nds) {
    auto n = ndpair.first;
    for (auto f : n->getAdjListOut()) {
      if (f->getFrom() != n) continue;
      if (_edgs.find(f) == _edgs.end()) {
        LOGTO(WARN, std::cerr) << "Edge " << f << " was not drawn, skipping...";
        continue;
      }

      double dTot = 0;
      for (const auto& seg : combEdgSegs[f])
        dTot += pathSegments[seg.first].geom.getLength();
      double tot = f->pl().getChilds().size();
      double step = dTot / tot;

      auto curSegIdx = 0;
      double curSegLen =
          pathSegments[combEdgSegs[f][curSegIdx].first].geom.getLength();
      double offset = 0;

      int i = 1;

      std::vector<LineNode*> nds{};
      auto prev = n->pl().getParent();
      for (auto e : f->pl().getChilds()) {
        if (e->getFrom() == prev)
          nds.push_back(e->getTo());
        else if (e->getTo() == prev)
          nds.push_back(e->getFrom());
        prev = nds.back();
      }

      nds.pop_back();

      for (auto nd : nds) {
        double progr = (step * i) / dTot;

        while (progr > offset + curSegLen) {
          curSegIdx++;
          offset += curSegLen;
          curSegLen =
              pathSegments[combEdgSegs[f][curSegIdx].first].geom.getLength();
        }

        double onSegProgr = progr - offset;
        assert(onSegProgr <= curSegLen);

        if (combEdgSegs[f][curSegIdx].second) {
          // using this segment in reverse direction
          pathSegments[combEdgSegs[f][curSegIdx].first].nodesOnSeg.push_back(
              {nd, 1 - onSegProgr});
        } else {
          pathSegments[combEdgSegs[f][curSegIdx].first].nodesOnSeg.push_back(
              {nd, onSegProgr});
        }

        i++;
      }
    }
  }

  std::vector<std::vector<LineEdge*>> segmentEdges(pathSegments.size());

  for (size_t segId = 0; segId < pathSegments.size(); segId++) {
    auto seg = pathSegments[segId];
    std::sort(seg.nodesOnSeg.begin(), seg.nodesOnSeg.end());

    std::set<shared::linegraph::LineEdge*> activeLineEdges;

    // first add new topological node at beginning of segment
    if (mm.find(seg.start) == mm.end())
      mm[seg.start] = target->addNd(seg.geom.front());

    auto from = mm[seg.start];
    double lastProgr = 0;

    for (size_t i = 0; i < seg.nodesOnSeg.size(); i++) {
      auto geom = seg.geom.getSegment(lastProgr, seg.nodesOnSeg[i].progression);

      auto to = seg.nodesOnSeg[i].n;

      if (m.find(to) == m.end()) {
        auto payload = to->pl();
        payload.setGeom(geom.getLine().back());
        m[to] = target->addNd(payload);
      }

      segmentEdges[segId].push_back(target->addEdg(from, m[to], geom));
      from = m[to];
      lastProgr = seg.nodesOnSeg[i].progression;
    }

    // add new topological node at end of segment
    if (mm.find(seg.end) == mm.end())
      mm[seg.end] = target->addNd(seg.geom.back());

    auto geom = seg.geom.getSegment(lastProgr, 1);
    segmentEdges[segId].push_back(target->addEdg(from, mm[seg.end], geom));
  }

  for (auto a : combEdgSegs) {
    auto combEdg = a.first;

    size_t childPointer = 0;
    auto prev = a.first->getFrom()->pl().getParent();

    for (auto segPairIt = a.second.rbegin(); segPairIt != a.second.rend();
         segPairIt++) {
      auto segPair = *segPairIt;
      size_t segId = segPair.first;
      bool reverse = segPair.second;

      auto sCopy = segmentEdges[segId];
      if (reverse) std::reverse(sCopy.begin(), sCopy.end());

      for (auto edge : sCopy) {
        auto from = reverse ? edge->getTo() : edge->getFrom();
        auto to = edge->getOtherNd(from);

        for (auto lo :
             combEdg->pl().getChilds()[childPointer]->pl().getLines()) {
          // TODO: direction
          edge->pl().addLine(lo.line, lo.direction);
        }

        auto curChldTo = combEdg->pl().getChilds()[childPointer]->getTo();
        if (curChldTo == prev)
          curChldTo = combEdg->pl().getChilds()[childPointer]->getFrom();

        if (to == m[curChldTo]) {
          prev = curChldTo;
          childPointer++;
        } else {
          for (auto lo :
               combEdg->pl().getChilds()[childPointer]->pl().getLines()) {
            to->pl().addLineNotServed(lo.line);
          }
        }
      }
    }
  }
}
// _____________________________________________________________________________
void Drawing::crumble() {
  _c = std::numeric_limits<double>::infinity();
  _nds.clear();
  _edgs.clear();
  _ndReachCosts.clear();
  _ndBndCosts.clear();
  _edgCosts.clear();
  _springCosts.clear();
}

// _____________________________________________________________________________
double Drawing::recalcBends(const CombNode* nd) {
  double c = 0;

  if (_nds.count(nd) == 0) return 0;
  auto gnd = _gg->getGrNdById(_nds.find(nd)->second);

  // TODO: implement this better

  for (auto e : nd->getAdjList()) {
    if (_edgs.count(e) == 0) {
      continue;  // dont count edge that havent been drawn
    }
    auto ge = _edgs.find(e)->second;

    size_t dirA = 0;
    for (; dirA < _gg->maxDeg(); dirA++) {
      if (gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.front())->getFrom() ||
          gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.front())->getTo() ||
          gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.back())->getFrom() ||
          gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.back())->getTo()) {
        break;
      }
    }

    assert(dirA < _gg->maxDeg());

    for (auto lo : e->pl().getChilds().front()->pl().getLines()) {
      for (auto f : nd->getAdjList()) {
        if (e == f) continue;
        if (_edgs.count(f) == 0) {
          continue;  // dont count edges that havent been drawn
        }
        auto gf = _edgs.find(f)->second;

        if (f->pl().getChilds().front()->pl().hasLine(lo.line)) {
          size_t dirB = 0;
          for (; dirB < _gg->maxDeg(); dirB++) {
            if (gnd->pl().getPort(dirB) ==
                    _gg->getGrEdgById(gf.front())->getFrom() ||
                gnd->pl().getPort(dirB) ==
                    _gg->getGrEdgById(gf.front())->getTo() ||
                gnd->pl().getPort(dirB) ==
                    _gg->getGrEdgById(gf.back())->getFrom() ||
                gnd->pl().getPort(dirB) ==
                    _gg->getGrEdgById(gf.back())->getTo()) {
              break;
            }
          }

          assert(dirB < _gg->maxDeg());

          c += _gg->getBendPen(dirA, dirB);
        }
      }
    }
  }

  return c / 2.0;
}

// _____________________________________________________________________________
bool Drawing::drawn(const CombEdge* ce) const { return _edgs.count(ce); }

// _____________________________________________________________________________
void Drawing::erase(CombEdge* ce) {
  _edgs.erase(ce);
  _c -= _edgCosts[ce];
  _edgCosts.erase(ce);

  _c -= _springCosts[ce];
  _springCosts.erase(ce);

  _c -= _ndBndCosts[ce->getFrom()];
  _c -= _ndBndCosts[ce->getTo()];

  // update bend costs
  _ndBndCosts[ce->getFrom()] = recalcBends(ce->getFrom());
  _ndBndCosts[ce->getTo()] = recalcBends(ce->getTo());

  _c += _ndBndCosts[ce->getTo()];
  _c += _ndBndCosts[ce->getFrom()];
}

// _____________________________________________________________________________
void Drawing::erase(CombNode* cn) {
  _nds.erase(cn);
  _c -= _ndReachCosts[cn];
  _c -= _ndBndCosts[cn];
  _ndReachCosts.erase(cn);
  _ndBndCosts.erase(cn);
}

// _____________________________________________________________________________
void Drawing::eraseFromGrid(const CombEdge* ce, BaseGraph* gg) {
  auto it = _edgs.find(ce);
  if (it == _edgs.end()) return;
  const auto& es = it->second;
  for (auto eid : es) {
    auto e = gg->getGrEdgById(eid);
    // TODO: remove const cast
    gg->unSettleEdg(const_cast<CombEdge*>(ce), e->getFrom()->pl().getParent(),
                    e->getTo()->pl().getParent());
  }
}

// _____________________________________________________________________________
void Drawing::applyToGrid(const CombEdge* ce, BaseGraph* gg) {
  auto it = _edgs.find(ce);
  if (it == _edgs.end()) return;
  const auto& es = it->second;

  for (auto eid : es) {
    auto e = gg->getGrEdgById(eid);
    gg->settleEdg(e->getFrom()->pl().getParent(), e->getTo()->pl().getParent(),
                  const_cast<CombEdge*>(ce), 0);
  }
}

// _____________________________________________________________________________
void Drawing::eraseFromGrid(const CombNode* nd, BaseGraph* gg) {
  gg->unSettleNd(const_cast<CombNode*>(nd));
}

// _____________________________________________________________________________
void Drawing::applyToGrid(const CombNode* nd, BaseGraph* gg) {
  gg->settleNd(const_cast<GridNode*>(gg->getGrNdById(_nds[nd])),
               const_cast<CombNode*>(nd));
}

// _____________________________________________________________________________
void Drawing::eraseFromGrid(BaseGraph* gg) {
  for (auto e : _edgs) eraseFromGrid(e.first, gg);
  for (auto nd : _nds) eraseFromGrid(nd.first, gg);
}

// _____________________________________________________________________________
void Drawing::applyToGrid(BaseGraph* gg) {
  for (auto nd : _nds) applyToGrid(nd.first, gg);
  for (auto e : _edgs) applyToGrid(e.first, gg);
}

// _____________________________________________________________________________
double Drawing::getEdgCost(const CombEdge* e) const {
  if (_edgCosts.count(e)) return _edgCosts.find(e)->second;
  return 0;
}

// _____________________________________________________________________________
double Drawing::getNdBndCost(const CombNode* n) const {
  if (_ndBndCosts.count(n)) return _ndBndCosts.find(n)->second;
  return 0;
}

// _____________________________________________________________________________
double Drawing::getNdReachCost(const CombNode* n) const {
  if (_ndReachCosts.count(n)) return _ndReachCosts.find(n)->second;
  return 0;
}
