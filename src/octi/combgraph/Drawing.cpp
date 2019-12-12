#include <iostream>
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/BezierCurve.h"
#include "util/graph/Dijkstra.h"
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;

using octi::combgraph::Drawing;
using octi::combgraph::Costs;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::CombEdge;
using util::graph::Dijkstra;
using util::geo::BezierCurve;
using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::gridgraph::GridNodePL;
using octi::gridgraph::GridEdgePL;

// _____________________________________________________________________________
double Drawing::score() const { return _c; }

// _____________________________________________________________________________
void Drawing::setGridGraph(const GridGraph* gg) {
  _gg = gg;
}

// _____________________________________________________________________________
Costs Drawing::fullScore() const {
  Costs ret{0, 0, 0, 0};

  for (auto c : _ndReachCosts) ret.move += c.second;
  for (auto c : _ndBndCosts) ret.bend += c.second;
  for (auto c : _edgCosts) ret.hop += c.second;
  for (auto c : _springCosts) ret.dense += c.second;

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
        assert(!ges[i + 1]->pl().isSecondary());
      }
    } else if (i == ges.size() - 1) {
      if (!_ndReachCosts.count(rev ? ce->getTo() : ce->getFrom())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[rev ? ce->getTo() : ce->getFrom()] = ge->pl().cost();
        _ndBndCosts[rev ? ce->getTo() : ce->getFrom()] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[rev ? ce->getTo() : ce->getFrom()] += ge->pl().cost();
        assert(!ges[i - 1]->pl().isSecondary());
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

  // variables named as in the paper
  int k = ce->pl().getChilds().size() - 1;

  double c = _gg->getPenalties().densityPen / (k);

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
PolyLine<double> Drawing::buildPolylineFromRes(
    const std::vector<std::pair<size_t, size_t>>& res) const {
  PolyLine<double> pl;
  for (auto revIt = res.rbegin(); revIt != res.rend(); revIt++) {
    auto f = _gg->getEdg(_gg->getGrNdById(revIt->first),
                         _gg->getGrNdById(revIt->second));
    // TODO check for isSecondary should not be needed, filtered out by draw()
    if (!f->pl().isSecondary()) {
      if (pl.getLine().size() > 0 &&
          dist(pl.getLine().back(), *f->getFrom()->pl().getGeom()) > 0) {
        BezierCurve<double> bc(pl.getLine().back(),
                               *f->getFrom()->pl().getParent()->pl().getGeom(),
                               *f->getFrom()->pl().getParent()->pl().getGeom(),
                               *f->getFrom()->pl().getGeom());

        for (auto p : bc.render(10).getLine()) pl << p;
      } else {
        pl << *f->getFrom()->pl().getParent()->pl().getGeom();
      }

      pl << *f->getFrom()->pl().getGeom();
      pl << *f->getTo()->pl().getGeom();
    }
  }

  if (res.size())
    pl << *_gg->getEdg(_gg->getGrNdById(res.front().first),
                       _gg->getGrNdById(res.front().second))
               ->getTo()
               ->pl()
               .getParent()
               ->pl()
               .getGeom();

  return pl;
}

// _____________________________________________________________________________
void Drawing::getLineGraph(LineGraph* target) const {
  std::map<LineNode*, LineNode*> m;

  for (auto ndpair : _nds) {
    auto n = ndpair.first;
    for (auto f : n->getAdjListOut()) {
      if (f->getFrom() != n) continue;
      auto poly = buildPolylineFromRes(_edgs.find(f)->second);
      double tot = f->pl().getChilds().size();
      double d = poly.getLength();
      double step = d / tot;

      int i = 0;

      auto pre = n->pl().getParent();

      for (auto e : f->pl().getChilds()) {
        auto from = e->getFrom();
        auto to = e->getTo();

        PolyLine<double> pl =
            poly.getSegment((step * i) / d, (step * (i + 1)) / d);

        if (from == pre) {
          pre = to;
        } else {
          pl.reverse();
          pre = from;
        }

        if (m.find(from) == m.end()) {
          shared::linegraph::LineNodePL payload = from->pl();
          payload.setGeom(pl.getLine().front());
          auto tfrom = target->addNd(payload);
          m[from] = tfrom;
        }

        if (m.find(to) == m.end()) {
          shared::linegraph::LineNodePL payload = to->pl();
          payload.setGeom(pl.getLine().back());
          auto tto = target->addNd(payload);
          m[to] = tto;
        }

        shared::linegraph::LineEdgePL payload = e->pl();
        payload.setPolyline(pl);
        target->addEdg(m[from], m[to], payload);

        i++;
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
  double c_0 = _gg->getPenalties().p_45 - _gg->getPenalties().p_135;
  double c_135 = _gg->getPenalties().p_45;
  double c_90 = _gg->getPenalties().p_45 - _gg->getPenalties().p_135 +
                _gg->getPenalties().p_90;
  double c_45 = c_0 + c_135;

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
    for (; dirA < 8; dirA++) {
      if (gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.front())->getFrom() ||
          gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.front())->getTo() ||
          gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.back())->getFrom() ||
          gnd->pl().getPort(dirA) == _gg->getGrEdgById(ge.back())->getTo()) {
        break;
      }
    }

    assert(dirA < 8);

    for (auto ro : e->pl().getChilds().front()->pl().getRoutes()) {
      for (auto f : nd->getAdjList()) {
        if (e == f) continue;
        if (_edgs.count(f) == 0) {
          continue;  // dont count edges that havent been drawn
        }
        auto gf = _edgs.find(f)->second;

        if (f->pl().getChilds().front()->pl().hasRoute(ro.route)) {
          size_t dirB = 0;
          for (; dirB < 8; dirB++) {
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

          assert(dirB < 8);

          int ang = (8 + (dirA - dirB)) % 8;
          if (ang > 4) ang = 8 - ang;
          ang = ang % 4;

          double mult = 1;

          // write corresponding cost to addC[j]
          if (ang == 0) c += mult * c_0;
          if (ang == 1) c += mult * c_45;
          if (ang == 2) c += mult * c_90;
          if (ang == 3) c += mult * c_135;
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
void Drawing::eraseFromGrid(const CombEdge* ce, GridGraph* gg) {
  auto it = _edgs.find(ce);
  if (it == _edgs.end()) return;
  const auto& es = it->second;
  for (auto eid : es) {
    auto e = gg->getGrEdgById(eid);
    gg->unSettleEdg(e->getFrom()->pl().getParent(),
                    e->getTo()->pl().getParent());
  }
}

// _____________________________________________________________________________
void Drawing::applyToGrid(const CombEdge* ce, GridGraph* gg) {
  auto it = _edgs.find(ce);
  if (it == _edgs.end()) return;
  const auto& es = it->second;

  for (auto eid : es) {
    auto e = gg->getGrEdgById(eid);
    gg->settleEdg(e->getFrom()->pl().getParent(), e->getTo()->pl().getParent(),
                  const_cast<CombEdge*>(ce));
  }
}

// _____________________________________________________________________________
void Drawing::eraseFromGrid(const CombNode* nd, GridGraph* gg) {
  gg->unSettleNd(const_cast<CombNode*>(nd));
}

// _____________________________________________________________________________
void Drawing::applyToGrid(const CombNode* nd, GridGraph* gg) {
  gg->settleNd(const_cast<GridNode*>(gg->getGrNdById(_nds[nd])),
               const_cast<CombNode*>(nd));
}

// _____________________________________________________________________________
void Drawing::eraseFromGrid(GridGraph* gg) {
  for (auto e : _edgs) eraseFromGrid(e.first, gg);
  for (auto nd : _nds) eraseFromGrid(nd.first, gg);
}

// _____________________________________________________________________________
void Drawing::applyToGrid(GridGraph* gg) {
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
