// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/combgraph/CombGraph.h"

using octi::combgraph::CombGraph;
using octi::transitgraph::TransitGraph;
using octi::transitgraph::TransitNode;

// _____________________________________________________________________________
CombGraph::CombGraph(const TransitGraph* g) {
  build(g);
  combineDeg2();
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
      addEdg(m[e->getFrom()], m[e->getTo()],
                     octi::combgraph::CombEdgePL(e));
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
