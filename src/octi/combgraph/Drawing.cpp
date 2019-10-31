#include <iostream>
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/BezierCurve.h"
#include "util/graph/Dijkstra.h"
using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;

using octi::combgraph::Drawing;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::CombEdge;
using util::graph::Dijkstra;
using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::gridgraph::GridNodePL;
using octi::gridgraph::GridEdgePL;

// _____________________________________________________________________________
double Drawing::score() const { return _c; }

// _____________________________________________________________________________
void Drawing::draw(CombEdge* ce, const GrEdgList& ges) {
  if (_edgs.count(ce)) _edgs[ce].clear();
  for (size_t i = 0; i < ges.size(); i++) {
    auto ge = ges[i];

      _nds[ce->getTo()] = ges.front()->getTo()->pl().getParent();
      _nds[ce->getFrom()] = ges.back()->getFrom()->pl().getParent();

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
      if (!_ndReachCosts.count(ce->getTo())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[ce->getTo()] = ge->pl().cost();
        _ndBndCosts[ce->getTo()] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[ce->getTo()] += ge->pl().cost();
      }
    } else if (i == ges.size() - 1) {
      if (!_ndReachCosts.count(ce->getFrom())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[ce->getFrom()] = ge->pl().cost();
        _ndBndCosts[ce->getFrom()] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[ce->getFrom()] += ge->pl().cost();
      }
    } else {
      _edgCosts[ce] += ge->pl().cost();
    }

      if (!ges[i]->pl().isSecondary()) _edgs[ce].push_back(ges[i]);
  }

  double sum = 0;
  for (auto a : _ndReachCosts) sum += a.second;
  for (auto a : _ndBndCosts) sum += a.second;
  for (auto a : _edgCosts) sum += a.second;

  if (fabs(sum - _c) > 0.0001) {
    // std::cerr << sum << " vs " << _c << std::endl;
    assert(false);
  }

  // for (auto a : _ndBndCosts) {
    // auto re = recalcBends(a.first);
    // if (fabs(re - a.second) > 0.0001) {
      // // std::cerr << "Recalc: " << re << " from edges: " << a.second << std::endl;
      // assert(false);
    // }
  // }
}

// _____________________________________________________________________________
const GridNode* Drawing::getGrNd(const CombNode* cn) { return _nds[cn]; }

// _____________________________________________________________________________
const std::vector<const GridEdge*>& Drawing::getGrEdgs(const CombEdge* ce) {
  return _edgs.find(ce)->second;
}
// _____________________________________________________________________________
PolyLine<double> Drawing::buildPolylineFromRes(
    const std::vector<const GridEdge*>& res) const {
  PolyLine<double> pl;
  for (auto revIt = res.rbegin(); revIt != res.rend(); revIt++) {
    auto f = *revIt;
    // TODO check for isSeconday should not be needed, filtered out by draw()
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

  if (res.size()) pl << *res.front()->getTo()->pl().getParent()->pl().getGeom();

  return pl;
}

// _____________________________________________________________________________
void Drawing::getTransitGraph(TransitGraph* target) const {
  std::map<TransitNode*, TransitNode*> m;

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
          auto payload = from->pl();
          payload.setGeom(pl.getLine().front());
          auto tfrom = target->addNd(payload);
          m[from] = tfrom;
        }

        if (m.find(to) == m.end()) {
          auto payload = to->pl();
          payload.setGeom(pl.getLine().back());
          auto tto = target->addNd(payload);
          m[to] = tto;
        }

        auto payload = e->pl();
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
}

// _____________________________________________________________________________
double Drawing::recalcBends(const CombNode* nd) {
  double c_0 = _gg->getPenalties().p_45 - _gg->getPenalties().p_135;
  double c_135 = _gg->getPenalties().p_45;
  double c_90 = _gg->getPenalties().p_45 - _gg->getPenalties().p_135 + _gg->getPenalties().p_90;
  double c_45 = c_0 + c_135;

  double c = 0;

  if (_nds.count(nd) == 0) return 0;
  auto gnd = _nds.find(nd)->second;

  // TODO: implement this better

  for (auto e : nd->getAdjList()) {
    if (_edgs.count(e) == 0) {
      continue;  // dont count edge that havent been drawn
    }
    auto ge = _edgs.find(e)->second;

    size_t dirA = 0;
    for (; dirA < 8; dirA++) {
      if (
          gnd->pl().getPort(dirA) == ge.front()->getFrom() ||
          gnd->pl().getPort(dirA) == ge.front()->getTo() ||
          gnd->pl().getPort(dirA) == ge.back()->getFrom() ||
          gnd->pl().getPort(dirA) == ge.back()->getTo()) {
        break;
      }
    }

    assert(dirA < 8);


    for (auto ro : e->pl().getChilds().front()->pl().getRoutes()) {
      for (auto f : nd->getAdjList()) {
        if (e == f) continue;
        if (_edgs.count(f) == 0) {
          continue;  // dont count edge that havent been drawn
        }
        auto gf = _edgs.find(f)->second;

        if (f->pl().getChilds().front()->pl().hasRoute(ro.route)) {
          size_t dirB = 0;
          for (; dirB < 8; dirB++) {
            if (
                gnd->pl().getPort(dirB) == gf.front()->getFrom() ||
                gnd->pl().getPort(dirB) == gf.front()->getTo() ||
                gnd->pl().getPort(dirB) == gf.back()->getFrom() ||
                gnd->pl().getPort(dirB) == gf.back()->getTo()) {
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

  return c / 2;
}

// _____________________________________________________________________________
bool Drawing::drawn(const CombEdge* ce) const {
  return _edgs.count(ce);
}

// _____________________________________________________________________________
void Drawing::erase(CombEdge* ce) {
  _edgs.erase(ce);
  _c -= _edgCosts[ce];
  _edgCosts.erase(ce);

  _c -= _ndBndCosts[ce->getFrom()];
  _c -= _ndBndCosts[ce->getTo()];

  // update bend costs
  _ndBndCosts[ce->getFrom()] = recalcBends(ce->getFrom());
  _ndBndCosts[ce->getTo()] = recalcBends(ce->getTo());

  _c += _ndBndCosts[ce->getTo()];
  _c += _ndBndCosts[ce->getFrom()];

  double sum = 0;
  for (auto a : _ndReachCosts) sum += a.second;
  for (auto a : _ndBndCosts) sum += a.second;
  for (auto a : _edgCosts) sum += a.second;

  if (fabs(sum - _c) > 0.0001) {
    std::cerr << sum << " vs " << _c << std::endl;
    assert(false);
  }

  assert(_edgs.count(ce) == 0);
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
  for (auto e : es) {
    gg->unSettleEdg(e->getFrom()->pl().getParent(),
                    e->getTo()->pl().getParent());
  }
}

// _____________________________________________________________________________
void Drawing::applyToGrid(const CombEdge* ce, GridGraph* gg) {
  auto it = _edgs.find(ce);
  if (it == _edgs.end()) return;
  const auto& es = it->second;

  for (auto e : es) {
    // if (e->pl().isSecondary()) continue;
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
  gg->settleNd(const_cast<GridNode*>(_nds[nd]), const_cast<CombNode*>(nd));
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
