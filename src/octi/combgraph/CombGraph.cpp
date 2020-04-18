// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/combgraph/CombGraph.h"

using octi::combgraph::CombGraph;
using octi::combgraph::EdgeOrdering;
using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;

// _____________________________________________________________________________
CombGraph::CombGraph(const LineGraph* g)
    : CombGraph(g, false) {}

// _____________________________________________________________________________
CombGraph::CombGraph(const LineGraph* g, bool collapse) {
  build(g);
  if (collapse) combineDeg2();
  writeEdgeOrdering();
}

// _____________________________________________________________________________
void CombGraph::build(const LineGraph* source) {
  auto nodes = source->getNds();

  std::map<LineNode*, CombNode*> m;

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
    size_t numLines = 0;
    for (auto e : n->getAdjList()) {
      numLines += e->pl().getChilds().front()->pl().getLines().size();
    }
    n->pl().setRouteNumber(numLines);
  }
}

// _____________________________________________________________________________
void CombGraph::combineDeg2() {
  for (auto n : *getNds()) {
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
void CombGraph::writeEdgeOrdering() {
  for (auto n : *getNds()) {
    n->pl().setEdgeOrdering(getEdgeOrderingForNode(n));
  }
}

// _____________________________________________________________________________
EdgeOrdering CombGraph::getEdgeOrderingForNode(CombNode* n) const {
  return getEdgeOrderingForNode(n, true);
}

// _____________________________________________________________________________
EdgeOrdering CombGraph::getEdgeOrderingForNode(CombNode* n,
    bool useOrigNextNode) const {
  EdgeOrdering order;
  for (auto e : n->getAdjList()) {
    auto r = e->pl().getChilds().front();
    util::geo::DPoint a = *n->pl().getGeom();

    util::geo::DPoint b;
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
      }
    }

    // get the angles
    double deg = util::geo::angBetween(a, b) - M_PI / 2;
    if (deg <= 0) deg += M_PI * 2;

    order.add(e, deg);
  }

  return order;
}
