// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/combgraph/CombGraph.h"

using octi::combgraph::CombGraph;
using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using octi::combgraph::EdgeOrdering;

// _____________________________________________________________________________
CombGraph::CombGraph(const TransitGraph* g) {
  build(g);
  combineDeg2();
  writeEdgeOrdering();
}

// _____________________________________________________________________________
void CombGraph::build(const TransitGraph* source) {
  auto nodes = source->getNds();

  std::map<TransitNode*, CombNode*> m;

  for (auto n : nodes) {
    CombNode* cn = addNd(n);
    m[n] = cn;
  }

  for (auto n : nodes) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      addEdg(m[e->getFrom()], m[e->getTo()], octi::combgraph::CombEdgePL(e));
    }
  }

  for (auto n : *getNds()) {
    size_t routes = 0;
    for (auto e : n->getAdjList()) {
      routes += e->pl().getChilds().front()->pl().getRoutes().size();
    }
    n->pl().setRouteNumber(routes);
  }
}

// _____________________________________________________________________________
void CombGraph::combineDeg2() {
  std::set<CombNode*> nodes = *getNds();
  for (auto n : nodes) {
    if (n->getAdjList().size() == 2) {
      CombEdge* a = n->getAdjList().front();
      CombEdge* b = n->getAdjList().back();

      // if this combination would turn our graph into a multigraph,
      // dont do it!
      if (getEdg(a->getOtherNd(n), b->getOtherNd(n))) continue;

      auto pl = a->pl();

      // a is the reference
      if (a->getTo() == n) {
        if (b->getTo() != n) {
          pl.getChilds().insert(pl.getChilds().end(),
                                b->pl().getChilds().begin(),
                                b->pl().getChilds().end());
        } else {
          pl.getChilds().insert(pl.getChilds().end(),
                                b->pl().getChilds().rbegin(),
                                b->pl().getChilds().rend());
        }

        pl.setPolyLine(PolyLine<double>(*a->getFrom()->pl().getGeom(),
                                        *b->getOtherNd(n)->pl().getGeom()));
        addEdg(a->getFrom(), b->getOtherNd(n), pl);
      } else {
        if (b->getTo() == n) {
          pl.getChilds().insert(pl.getChilds().begin(),
                                b->pl().getChilds().begin(),
                                b->pl().getChilds().end());
        } else {
          pl.getChilds().insert(pl.getChilds().begin(),
                                b->pl().getChilds().rbegin(),
                                b->pl().getChilds().rend());
        }

        pl.setPolyLine(PolyLine<double>(*b->getOtherNd(n)->pl().getGeom(),
                                        *a->getTo()->pl().getGeom()));
        addEdg(b->getOtherNd(n), a->getTo(), pl);
      }

      delNd(n);
    }
  }
}

// _____________________________________________________________________________
void CombGraph::getTransitGraph(TransitGraph* target) const {
  std::map<TransitNode*, TransitNode*> m;

  for (auto n : getNds()) {
    for (auto f : n->getAdjListOut()) {
      if (f->getFrom() != n) continue;
      if (f->pl().getGeneration() < 0) continue;
      double tot = f->pl().getChilds().size();
      double d = f->pl().getPolyLine().getLength();
      double step = d / tot;

      int i = 0;

      auto pre = n->pl().getParent();

      for (auto e : f->pl().getChilds()) {
        auto from = e->getFrom();
        auto to = e->getTo();

        PolyLine<double> pl = f->pl().getPolyLine().getSegment(
            (step * i) / d, (step * (i + 1)) / d);

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
        payload.setGeneration(f->pl().getGeneration());
        target->addEdg(m[from], m[to], payload);

        i++;
      }
    }
  }
}

// _____________________________________________________________________________
void CombGraph::writeEdgeOrdering() {
  for (auto n : *getNds()) {
    n->pl().setEdgeOrdering(getEdgeOrderingForNode(n));
  }
}

// _____________________________________________________________________________
EdgeOrdering CombGraph::getEdgeOrderingForNode(CombNode* n) const {
  return getEdgeOrderingForNode(n, true, std::map<CombNode*, DPoint>());
}

// _____________________________________________________________________________
EdgeOrdering CombGraph::getEdgeOrderingForNode(
    CombNode* n, bool useOrigNextNode,
    const std::map<CombNode*, DPoint>& newPos) const {
  EdgeOrdering order;
  for (auto e : n->getAdjList()) {
    auto r = e->pl().getChilds().front();
    DPoint a = *n->pl().getGeom();
    if (newPos.find(n) != newPos.end()) a = newPos.find(n)->second;

    DPoint b;
    if (useOrigNextNode) {
      b = *r->getOtherNd(n->pl().getParent())->pl().getGeom();
    } else {
      auto other = e->getOtherNd(n);
      if (e->pl().getGeom()->size() > 2) {
        if (e->getTo() == n) {
          b = e->pl().getGeom()->at(e->pl().getGeom()->size() - 2);
        } else {
          b = e->pl().getGeom()->at(1);
        }
      } else {
        b = *other->pl().getGeom();
        if (newPos.find(other) != newPos.end()) b = newPos.find(other)->second;
      }
    }

    // get the angles
    double deg = util::geo::angBetween(a, b);

    order.add(e, deg);
  }

  return order;
}

// _____________________________________________________________________________
size_t CombGraph::changesTopology(
    CombNode* nOrig, DPoint p,
    const std::map<CombNode*, DPoint>& newPos) const {
  // collect the affected nodes
  size_t ret = 0;
  std::set<CombNode*> aff;
  auto newPosA = newPos;
  newPosA[nOrig] = p;
  for (auto e : nOrig->getAdjList()) {
    if (newPos.find(e->getFrom()) != newPos.end()) aff.insert(e->getFrom());
    if (newPos.find(e->getTo()) != newPos.end()) aff.insert(e->getTo());
  }

  for (auto n : aff) {
    if (!getEdgeOrderingForNode(n, false, newPosA)
             .equals(getEdgeOrderingForNode(n))) {
      ret += 1;
    }
  }

  return ret;
}
