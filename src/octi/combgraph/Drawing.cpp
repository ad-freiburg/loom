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
using octi::combgraph::GrPath;
using octi::combgraph::Score;
using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using util::geo::BezierCurve;
using util::graph::Dijkstra;

// _____________________________________________________________________________
double Drawing::score() const {
  return _c + violations() * basegraph::SOFT_INF;
}

// _____________________________________________________________________________
double Drawing::rawScore() const {
  return _c + violations() * basegraph::SOFT_INF;
}

// _____________________________________________________________________________
uint64_t Drawing::violations() const { return _violations; }

// _____________________________________________________________________________
void Drawing::setBaseGraph(const BaseGraph* gg) { _gg = gg; }

// _____________________________________________________________________________
Score Drawing::fullScore() const {
  Score ret{0, 0, 0, 0, 0, 0};

  for (auto c : _ndReachCosts) ret.move += c.second;
  for (auto c : _ndBndCosts) ret.bend += c.second;
  for (auto c : _edgCosts) ret.hop += c.second;
  for (auto c : _springCosts) ret.dense += c.second;
  ret.full = _c + basegraph::SOFT_INF * violations();
  ret.violations = violations();

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

    double edgeCost = ge->pl().cost();
    if (edgeCost >= basegraph::SOFT_INF) {
      int vios = edgeCost / basegraph::SOFT_INF;
      edgeCost -= vios * basegraph::SOFT_INF;
      _vios[ce]++;
      _violations++;
    }

    _c += edgeCost;

    if (i == 0) {
      if (!_ndReachCosts.count(rev ? ce->getFrom() : ce->getTo())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[rev ? ce->getFrom() : ce->getTo()] = edgeCost;
        _ndBndCosts[rev ? ce->getFrom() : ce->getTo()] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[rev ? ce->getFrom() : ce->getTo()] += edgeCost;
      }
    } else if (i == ges.size() - 1) {
      if (!_ndReachCosts.count(rev ? ce->getTo() : ce->getFrom())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[rev ? ce->getTo() : ce->getFrom()] = edgeCost;
        _ndBndCosts[rev ? ce->getTo() : ce->getFrom()] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[rev ? ce->getTo() : ce->getFrom()] += edgeCost;
      }
    } else {
      if (!ge->pl().isSecondary()) l++;
      _edgCosts[ce] += edgeCost;
    }

    if (rev) {
      auto e = _gg->getEdg(ges[ges.size() - 1 - i]->getTo(),
                           ges[ges.size() - 1 - i]->getFrom());

      if (!e->pl().isSecondary()) {
        _edgs[ce].push_back(
            {e->getFrom()->pl().getId(), e->getTo()->pl().getId()});
      }
    } else {
      if (!ges[i]->pl().isSecondary()) {
        _edgs[ce].push_back(
            {ges[i]->getFrom()->pl().getId(), ges[i]->getTo()->pl().getId()});

        assert(_gg->getEdg(ges[i]->getFrom(), ges[i]->getTo()) == ges[i]);

        assert(
            _gg->getGrNdById(ges[i]->getFrom()->pl().getId())->pl().getId() ==
            ges[i]->getFrom()->pl().getId());

        assert(_gg->getGrNdById(ges[i]->getFrom()->pl().getId()) ==
               ges[i]->getFrom());

        assert(_gg->getGrEdgById({ges[i]->getFrom()->pl().getId(),
                                  ges[i]->getTo()->pl().getId()}) == ges[i]);
      }
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
  // This method extracts a line graph from the drawing, written to "target"
  // The challenge here is that because of the constraint relaxation, it may be
  // that image paths overlap. In this case, we have to merge overlapping
  // segments into a new segment in the target linegraph that wasn't originally
  // present in the input linegraph

  std::map<LineNode*, LineNode*> m;
  std::map<GridNode*, LineNode*> mm;

  std::vector<Segment> pathSegs;
  std::map<const CombEdge*, std::vector<std::pair<size_t, bool>>> cEdgSeg;

  std::map<const GridEdge*, size_t> grEdgToPathSeg;

  // settle grid nodes, _nds contains a mapping of input comb edges to
  // grid node ids
  for (auto ndpair : _nds) {
    auto combNd = ndpair.first;
    for (auto f : combNd->getAdjListOut()) {
      // go over each adjacent edge's image path and add nodes to the target
      // graph for the image's end and start node

      if (f->getFrom() != combNd) continue;
      if (_edgs.find(f) == _edgs.end()) {
        LOGTO(WARN, std::cerr) << "Edge " << f << " was not drawn, skipping...";
        continue;
      }

      // the image path...
      auto pth = _edgs.find(f)->second;
      assert(_gg->getGrEdgById(pth.back()));
      assert(_gg->getGrEdgById(pth.front()));
      // ... and it's from and to grid nodes. We can be sure that that
      // the path indeed goes in the direction of the comb edge, as the path is
      // explicitely reverted during the drawing if that was not the case
      GridNode* fr = _gg->getGrEdgById(pth.back())->getFrom()->pl().getParent();
      GridNode* to = _gg->getGrEdgById(pth.front())->getTo()->pl().getParent();

      // if we haven't seen the grid nodes...
      if (mm.find(fr) == mm.end()) {
        // ... retrieve a copy the corresponding PL of the comb node ...
        auto pl = combNd->pl().getParent()->pl();
        // ... update the geometry to the grid node geometry ...
        pl.setGeom(*fr->pl().getGeom());
        // ... and update the lookup maps
        mm[fr] = target->addNd(pl);
        m[combNd->pl().getParent()] = mm[fr];
      }

      // (as above)
      if (mm.find(to) == mm.end()) {
        auto pl = f->getTo()->pl().getParent()->pl();
        pl.setGeom(*to->pl().getGeom());
        mm[to] = target->addNd(pl);
        m[f->getTo()->pl().getParent()] = mm[to];
      }
    }
  }

  // build segments per path
  for (auto ndpair : _nds) {
    auto n = ndpair.first;
    for (auto f : n->getAdjListOut()) {
      if (f->getFrom() != n) continue;
      if (_edgs.find(f) == _edgs.end()) continue; // edge was not drawn

      auto path = _edgs.find(f)->second;

      std::set<CombEdge*> curResEdgs;

      for (size_t i = 0; i < path.size(); i++) {
        // the current edge id of the path
        auto pathEdgId = path[i];
        const GridEdge* grEdg = _gg->getGrEdgById(pathEdgId);
        assert(grEdg);
        assert(_gg->getResEdgsDirInd(grEdg).size());

        // the same grid edge in the other direction
        const GridEdge* grCounterEdg =
            _gg->getGrEdgById({pathEdgId.second, pathEdgId.first});

        // if we have already encountered this grid edge before, it hosts
        // multiple comb edges
        if (grEdgToPathSeg.count(grEdg)) {
          size_t exiSeg = grEdgToPathSeg[grEdg];
          if (cEdgSeg[f].size() == 0 || cEdgSeg[f].back().first != exiSeg) {
            // add this segment to the path
            cEdgSeg[f].push_back({exiSeg, false});
            curResEdgs = _gg->getResEdgsDirInd(grEdg);
            continue;
          }
        } else if (grEdgToPathSeg.count(grCounterEdg)) {
          // as above
          size_t exiSeg = grEdgToPathSeg[grCounterEdg];
          if (cEdgSeg[f].size() == 0 || cEdgSeg[f].back().first != exiSeg) {
            cEdgSeg[f].push_back({exiSeg, true});
            curResEdgs = _gg->getResEdgsDirInd(grCounterEdg);
            continue;
          }
        }

        // if we have not encountered it yet, check the resident edges on this
        // segment
        auto newResEdgs = _gg->getResEdgsDirInd(grEdg);
        assert(newResEdgs.size());

        if (newResEdgs == curResEdgs) {
          // if it is the same as before, extend the current segment
          pathSegs[cEdgSeg[f].back().first].start =
              grEdg->getFrom()->pl().getParent();
          pathSegs[cEdgSeg[f].back().first].path.push_back(pathEdgId);
        } else {
          // if not, start a new segment
          pathSegs.push_back({grEdg->getFrom()->pl().getParent(),
                              grEdg->getTo()->pl().getParent(),
                              {pathEdgId},
                              {},
                              {},
                              newResEdgs});
          cEdgSeg[f].push_back({pathSegs.size() - 1, false});
          curResEdgs = newResEdgs;
        }
        grEdgToPathSeg[grEdg] = cEdgSeg[f].back().first;
      }
    }
  }

  // write geometries of segments
  for (auto& segment : pathSegs) segment.geom = _gg->geomFromPath(segment.path);

  // add nodes to segments
  for (auto ndpair : _nds) {
    auto n = ndpair.first;
    for (auto f : n->getAdjListOut()) {
      if (f->getFrom() != n) continue;
      if (_edgs.find(f) == _edgs.end()) continue;

      double dTot = 0;
      for (const auto& seg : cEdgSeg[f])
        dTot += pathSegs[seg.first].geom.getLength();
      double tot = f->pl().getChilds().size();
      double step = dTot / tot;

      int i = 1;

      std::vector<LineNode*> nds;
      auto prev = n->pl().getParent();
      for (auto e : f->pl().getChilds()) {
        if (e->getFrom() == prev)
          nds.push_back(e->getTo());
        else if (e->getTo() == prev)
          nds.push_back(e->getFrom());
        prev = nds.back();
      }

      nds.pop_back();

      auto curSegIdx = cEdgSeg[f].size() - 1;
      double curSegLen = pathSegs[cEdgSeg[f][curSegIdx].first].geom.getLength();
      double offset = 0;

      for (auto nd : nds) {
        double progr = step * i;

        while (progr > offset + curSegLen) {
          curSegIdx--;
          offset += curSegLen;
          curSegLen = pathSegs[cEdgSeg[f][curSegIdx].first].geom.getLength();
        }

        double onSegProgr = (progr - offset) / curSegLen;

        if (cEdgSeg[f][curSegIdx].second) {
          // using this segment in reverse direction
          pathSegs[cEdgSeg[f][curSegIdx].first].nodes.push_back(
              {nd, 1 - onSegProgr});
        } else {
          pathSegs[cEdgSeg[f][curSegIdx].first].nodes.push_back(
              {nd, onSegProgr});
        }

        i++;
      }
    }
  }

  std::vector<std::vector<LineEdge*>> segmentEdges(pathSegs.size());

  for (size_t segId = 0; segId < pathSegs.size(); segId++) {
    auto seg = pathSegs[segId];
    std::sort(seg.nodes.begin(), seg.nodes.end());

    std::set<shared::linegraph::LineEdge*> activeLineEdges;

    // first add new topological node at beginning of segment
    if (mm.find(seg.start) == mm.end()) {
      mm[seg.start] = target->addNd(seg.geom.front());
    }

    auto from = mm[seg.start];
    double lastProgr = 0;

    for (size_t i = 0; i < seg.nodes.size(); i++) {
      auto geom = seg.geom.getSegment(lastProgr, seg.nodes[i].progr);

      auto to = seg.nodes[i].n;

      if (m.find(to) == m.end()) {
        auto payload = to->pl();
        payload.setGeom(geom.getLine().back());
        m[to] = target->addNd(payload);
      }

      segmentEdges[segId].push_back(target->addEdg(from, m[to], geom));
      from = m[to];
      lastProgr = seg.nodes[i].progr;
    }

    // add new topological node at end of segment
    if (mm.find(seg.end) == mm.end())
      mm[seg.end] = target->addNd(seg.geom.back());

    auto geom = seg.geom.getSegment(lastProgr, 1);
    segmentEdges[segId].push_back(target->addEdg(from, mm[seg.end], geom));
  }

  // now go over all original comb edge segments and add line informations
  // to the new line graph edges
  for (const auto& a : cEdgSeg) {
    auto combEdg = a.first;

    size_t childPtr = 0;
    auto prev = a.first->getFrom()->pl().getParent();

    // reverse iterator because paths are always given in reverse order
    for (auto it = a.second.rbegin(); it != a.second.rend(); it++) {
      auto segPair = *it;
      size_t segId = segPair.first;
      bool reverse = segPair.second;

      auto sCopy = segmentEdges[segId];
      if (reverse) std::reverse(sCopy.begin(), sCopy.end());

      for (auto edge : sCopy) {
        auto from = reverse ? edge->getTo() : edge->getFrom();
        auto to = edge->getOtherNd(from);

        assert(childPtr < combEdg->pl().getChilds().size());
        auto child = combEdg->pl().getChilds()[childPtr];

        edge->pl().setComponent(child->pl().getComponent());
        edge->getFrom()->pl().setComponent(child->pl().getComponent());
        edge->getTo()->pl().setComponent(child->pl().getComponent());

        for (auto lo : child->pl().getLines()) {
          // directions will be handled below by nodeRpl
          edge->pl().addLine(lo.line, lo.direction);
        }

        assert(child->getTo() == prev || child->getFrom() == prev);

        auto curChldTo = child->getTo();
        if (curChldTo == prev) curChldTo = child->getFrom();
        auto curChldFr = child->getOtherNd(curChldTo);

        if (from == m[curChldFr]) {
          // if this edge is the beginning of an original child edge, replace
          // it in the node. All other nodes will not be touched (except for
          // line not served info) by this original comb edge
          LineGraph::edgeRpl(m[curChldFr], child, edge);
        }

        if (to == m[curChldTo]) {
          // if this edge is the end of an original child edge, replace
          // it in the node. All other nodes will not be touched (except for
          // lines not served info) by this original comb edge
          LineGraph::edgeRpl(m[curChldTo], child, edge);
          prev = curChldTo;

          // consider the next child in the original comb edge
          childPtr++;
        } else {
          for (auto lo : child->pl().getLines()) {
            to->pl().addLineNotServed(lo.line);
          }
        }

        // node replacement in the altered edge, where directions are now given
        // w.r.t. to the original end node of the child edge. These are replaced
        // by the respective nodes of the new edge
        LineGraph::nodeRpl(edge, curChldFr, from);
        LineGraph::nodeRpl(edge, curChldTo, to);
      }
    }
  }
}
// _____________________________________________________________________________
void Drawing::crumble() {
  _c = std::numeric_limits<double>::infinity();
  _violations = 0;
  _nds.clear();
  _edgs.clear();
  _ndReachCosts.clear();
  _ndBndCosts.clear();
  _edgCosts.clear();
  _vios.clear();
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
const std::map<const CombEdge*, GrPath>& Drawing::getEdgPaths() const {
  return _edgs;
}

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

  _violations -= _vios[ce];
  _vios.erase(ce);
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
