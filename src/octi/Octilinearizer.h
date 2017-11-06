// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_OCTILINEARIZER_H_
#define OCTI_OCTILINEARIZER_H_

#include <unordered_set>
#include "octi/graph/Graph.h"
#include "octi/gridgraph/GridGraph.h"

namespace octi {

typedef util::graph::Graph<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombGraph;
typedef util::graph::Node<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombNode;
typedef util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombEdge;

using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::gridgraph::Penalties;
typedef octi::graph::Graph TransitGraph;
typedef util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL> TransitNode;
typedef util::graph::Edge<octi::graph::NodePL, octi::graph::EdgePL> TransitEdge;

struct GridHeur {
  GridHeur(octi::gridgraph::GridGraph* g, GridNode* from, const std::unordered_set<GridNode*>& to) : g(g), from(from), to(0) {
    if (to.size() == 1) {
      this->to = *to.begin();
    }

    for (auto n : to) {
      auto coords = g->getNodeCoords(n);

      for (size_t i = 0; i < 8; i++) {
        auto neigh = g->getNeighbor(coords.first, coords.second, i);
        if (neigh && to.find(neigh) == to.end()) {
          hull.push_back(g->getNodeCoords(n));
          break;
        }
      }
    }

    std::cerr << "to size: " << to.size() << ", hull size: " << hull.size() << std::endl;
  }

  double operator()(GridNode* from, const std::unordered_set<GridNode*>& to) {
    //return 0;

    if (to.find(from->pl().getParent()) != to.end()) return 0;

    size_t ret = std::numeric_limits<size_t>::max();
    auto xy = g->getNodeCoords(from->pl().getParent());

    for (auto t : hull) {
      size_t temp = g->heurCost(xy.first, xy.second, t.first, t.second);
      if (temp < ret) ret = temp;
    }

    return ret;
  }

  octi::gridgraph::GridGraph* g;
  GridNode* from;
  GridNode* to;
  std::vector<std::pair<size_t, size_t> > hull;
};

class Octilinearizer {
 public:
  Octilinearizer() {}

  TransitGraph draw(TransitGraph* g, const Penalties& pens);
 private:
  void normalizeCostVector(double* vec) const;
  double getMaxDis(CombNode* to, CombEdge* e);
  void removeEdgesShorterThan(TransitGraph* g, double d);
  void combineDeg2(CombGraph* g);
  void buildTransitGraph(CombGraph* source, TransitGraph* target);
  void writeEdgeOrdering(CombGraph* target);
  void buildCombGraph(TransitGraph* source, CombGraph* target);
  void rotate(CombGraph* g, Point center, double deg);
  void buildPolylineFromRes(const std::list<GridEdge*>& l, PolyLine& res);
  double getCostFromRes(const std::list<GridEdge*>& l);
  double addResidentEdges(gridgraph::GridGraph* g, CombEdge* e, const std::list<GridEdge*>& res);
};

}  // namespace octt

#endif  // OCTI_OCTILINEARIZER_H_
