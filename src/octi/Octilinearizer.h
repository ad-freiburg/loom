// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_OCTILINEARIZER_H_
#define OCTI_OCTILINEARIZER_H_

#include <vector>
#include <unordered_set>
#include "octi/graph/Graph.h"
#include "util/graph/Dijkstra.h"
#include "octi/gridgraph/GridGraph.h"

namespace octi {

typedef util::graph::UndirGraph<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombGraph;
typedef util::graph::Node<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombNode;
typedef util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombEdge;

using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::gridgraph::Penalties;
typedef octi::graph::Graph TransitGraph;
typedef util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL> TransitNode;
typedef util::graph::Edge<octi::graph::NodePL, octi::graph::EdgePL> TransitEdge;

struct GraphMeasures {
  double maxNodeDist;
  double minNodeDist;
  double avgNodeDist;
  double maxEdgeLength;
  double minEdgeLength;
  double avgEdgeLength;
  size_t maxDeg;
};

struct GridCost : public util::graph::Dijkstra::CostFunc<gridgraph::NodePL, gridgraph::EdgePL, double> {
  virtual double operator()(const GridNode* from, const GridEdge* e,
                       const GridNode* to) const {
    UNUSED(from);
    UNUSED(to);
    return e->pl().cost();
  }

  virtual double inf() const { return std::numeric_limits<double>::infinity(); }
};

struct GridHeur : public util::graph::Dijkstra::HeurFunc<gridgraph::NodePL, gridgraph::EdgePL, double> {
  GridHeur(octi::gridgraph::GridGraph* g, GridNode* from,
           const std::set<GridNode*>& to)
      : g(g), from(from), to(0) {
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
  }

  double operator()(const GridNode* from, const std::set<GridNode*>& to) const {
    // return 0;

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

  TransitGraph draw(TransitGraph* g, GridGraph** gg, const Penalties& pens);

 private:
  void normalizeCostVector(double* vec) const;
  double getMaxDis(CombNode* to, CombEdge* e, double gridSize);
  void removeEdgesShorterThan(TransitGraph* g, double d);
  void combineDeg2(CombGraph* g);
  void buildTransitGraph(CombGraph* source, TransitGraph* target);
  void writeEdgeOrdering(CombGraph* target);
  graph::EdgeOrdering getEdgeOrderingForNode(CombNode* n) const;
  graph::EdgeOrdering getEdgeOrderingForNode(
      CombNode* n, bool useOrigNextNode,
      const std::map<CombNode*, DPoint>& newPos) const;
  void buildCombGraph(TransitGraph* source, CombGraph* target);
  void buildPolylineFromRes(const std::vector<GridEdge*>& l,
                            PolyLine<double>& res);
  double getCostFromRes(const std::vector<GridEdge*>& l);
  void addResidentEdges(gridgraph::GridGraph* g, CombEdge* e,
                          const std::vector<GridEdge*>& res);
  size_t changesTopology(CombNode* n, DPoint p,
                         const std::map<CombNode*, DPoint>& newPos) const;
};

}  // namespace octt

#endif  // OCTI_OCTILINEARIZER_H_
