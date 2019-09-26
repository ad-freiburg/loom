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
void Drawing::draw(CombEdge* ce, const GrEdgList& ges, bool rev) {
  std::cerr << "Drawing " << ce << std::endl;
  _edgs[ce].clear();
  for (size_t i = 0; i < ges.size(); i++) {
    auto ge = ges[i];

    if (rev) {
      _nds[ce->getFrom()] = ges.front()->getFrom()->pl().getParent();
      _nds[ce->getTo()] =  ges.back()->getTo()->pl().getParent();
    } else {
      _nds[ce->getTo()] = ges.front()->getFrom()->pl().getParent();
      _nds[ce->getFrom()] =  ges.back()->getTo()->pl().getParent();
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
      std::cerr << "0: " << ge->pl().cost() << std::endl;
      if (!_ndReachCosts.count(rev ? ce->getTo() : ce->getFrom())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[rev ? ce->getTo() : ce->getFrom()] = ge->pl().cost();
        _ndBndCosts[{rev ? ce->getTo() : ce->getFrom(), ce}] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[{rev ? ce->getTo() : ce->getFrom(), ce}] += ge->pl().cost();
      }
    } else if (i == ges.size() - 1) {
      std::cerr << "last: " << ge->pl().cost() << std::endl;
      if (!_ndReachCosts.count(rev ? ce->getFrom() : ce->getTo())) {
        // if the node was not settled before, this is the node move cost
        _ndReachCosts[rev ? ce->getFrom() : ce->getTo()] = ge->pl().cost();
        _ndBndCosts[{rev ? ce->getFrom() : ce->getTo(), ce}] = 0;
      } else {
        // otherwise it is the reach cost belonging to the edge
        _ndBndCosts[{rev ? ce->getFrom() : ce->getTo(), ce}] += ge->pl().cost();
      }
    } else {
      std::cerr << "in: " << ge->pl().cost() << std::endl;
        _edgCosts[ce] += ge->pl().cost();
    }

    if (!ge->pl().isSecondary()) {
      if (rev) {
        _edgs[ce].push_back(_gg->getEdg(ges[ges.size() - 1 - i]->getTo(),
                                        ges[ges.size() - 1 - i]->getFrom()));
      } else {
        _edgs[ce].push_back(ges[i]);
      }
    }
  }

  double sum = 0;
  for (auto a : _ndReachCosts) sum += a.second;
  for (auto a : _ndBndCosts) sum += a.second;
  for (auto a : _edgCosts) sum += a.second;

  if (fabs(sum - _c) > 0.0001) {
    std::cerr << sum << " vs " << _c << std::endl;
    assert(false);
  }
}

// _____________________________________________________________________________
const GridNode* Drawing::getGrNd(const CombNode* cn) {
  return _nds[cn];
}

// _____________________________________________________________________________
const std::vector<const GridEdge*>& Drawing::getGrEdgs(const CombEdge* ce) {
  return _edgs[ce];
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
void Drawing::recalcBends(CombNode* nd) {
  std::cerr << "Recalcing bend costs at node " << nd << std::endl;
}

// _____________________________________________________________________________
void Drawing::erase(CombEdge* ce) {
  _edgs.erase(ce);
  _c -= _edgCosts[ce];
  _edgCosts.erase(ce);

  // update bend costs
  recalcBends(ce->getFrom());
  recalcBends(ce->getTo());
}

// _____________________________________________________________________________
void Drawing::erase(CombNode* cn) {
  _nds.erase(cn);
  _c -= _ndReachCosts[cn];
  // _c -= _ndBndCosts[cn];
  _ndReachCosts.erase(cn);
  // _ndBndCosts.erase(cn);
}
