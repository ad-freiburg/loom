// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_OCTILINEARIZER_H_
#define OCTI_OCTILINEARIZER_H_

#include <unordered_map>
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
  GridHeur(octi::gridgraph::GridGraph* g, GridNode* from, const std::unordered_map<GridNode*, bool>& to) : g(g), from(from), to(0) {
    if (to.size() == 1) {
      this->to = to.begin()->first;
    }
  /*
   *     minX = std::numeric_limits<size_t>::max();
   *     minY = std::numeric_limits<size_t>::max();
   *     maxX = 0;
   *     maxY = 0;
   *
   *     for (auto n : to) {
   *       auto xy = g->getNodeCoords(n.first->pl().getParent());
   *       if (xy.first < minX) minX = xy.first;
   *       if (xy.second < minY) minY = xy.second;
   *       if (xy.first > maxX) maxX = xy.first;
   *       if (xy.second > maxY) maxY = xy.second;
   *     }
   *
   *     for (auto n : to) {
   *       auto xy = g->getNodeCoords(n.first->pl().getParent());
   *       if (xy.first == minX || xy.first == maxX || xy.second == minY || xy.second == maxY) hull[n.first] = true;
   *     }
   */
  }

  double operator()(GridNode* from, const std::unordered_map<GridNode*, bool>& to) {
    size_t ret = std::numeric_limits<size_t>::max();
    auto xy = g->getNodeCoords(from->pl().getParent());

    for (auto t : to) {
      auto txy = g->getNodeCoords(t.first);
      size_t temp = g->heurCost(xy.first, xy.second, txy.first, txy.second);
      if (temp < ret) ret = temp;
    }

    return ret;
  }

  octi::gridgraph::GridGraph* g;
  GridNode* from;
  GridNode* to;
  std::unordered_map<GridNode*, bool> hull;
  size_t minX;
  size_t minY;
  size_t maxX;
  size_t maxY;
};

class Octilinearizer {
 public:
  Octilinearizer() {}

  TransitGraph draw(TransitGraph* g, const Penalties& pens);
 private:
  double getMaxDis(CombNode* to, CombEdge* e);
  void removeEdgesShorterThan(TransitGraph* g, double d);
  void combineDeg2(CombGraph* g);
  void buildTransitGraph(CombGraph* source, TransitGraph* target);
  void writeEdgeOrdering(CombGraph* target);
  void buildCombGraph(TransitGraph* source, CombGraph* target);
  void rotate(CombGraph* g, Point center, double deg);
};

}  // namespace octt


#endif  // OCTI_OCTILINEARIZER_H_
