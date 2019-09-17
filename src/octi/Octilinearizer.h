// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_OCTILINEARIZER_H_
#define OCTI_OCTILINEARIZER_H_

#include <unordered_set>
#include <vector>
#include "octi/combgraph/CombGraph.h"
#include "octi/gridgraph/GridGraph.h"
#include "shared/transitgraph/TransitGraph.h"
#include "util/graph/Dijkstra.h"

namespace octi {

using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::gridgraph::GridNodePL;
using octi::gridgraph::GridEdgePL;
using octi::gridgraph::Penalties;
using octi::gridgraph::NodeCost;
using octi::gridgraph::GrEdgList;

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdge;

using octi::combgraph::EdgeOrdering;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::CombEdge;

using util::graph::Dijkstra;

using octi::gridgraph::GrEdgList;

// comparator for nodes, based on degree
struct NodeCmp {
  bool operator()(CombNode* a, CombNode* b) {
    // return a->getAdjList().size() < b->getAdjList().size();
    return a->getAdjList().size() < b->getAdjList().size() ||
           (a->getAdjList().size() == b->getAdjList().size() &&
            a->pl().getRouteNumber() < b->pl().getRouteNumber());
  }
};

typedef std::priority_queue<CombNode*, std::vector<CombNode*>, NodeCmp> NodePQ;

struct GraphMeasures {
  double maxNodeDist;
  double minNodeDist;
  double avgNodeDist;
  double maxEdgeLength;
  double minEdgeLength;
  double avgEdgeLength;
  size_t maxDeg;
};

struct GridCost : public Dijkstra::CostFunc<GridNodePL, GridEdgePL, double> {
  virtual double operator()(const GridNode* from, const GridEdge* e,
                            const GridNode* to) const {
    UNUSED(from);
    UNUSED(to);
    return e->pl().cost();
  }

  virtual double inf() const { return std::numeric_limits<double>::infinity(); }
};

struct GridHeur : public Dijkstra::HeurFunc<GridNodePL, GridEdgePL, double> {
  GridHeur(GridGraph* g, GridNode* from, const std::set<GridNode*>& to)
      : g(g), from(from), to(0) {
    if (to.size() == 1) {
      this->to = *to.begin();
    }

    for (auto n : to) {
      for (size_t i = 0; i < 8; i++) {
        auto neigh = g->getNeighbor(n->pl().getX(), n->pl().getY(), i);
        if (neigh && to.find(neigh) == to.end()) {
          hull.push_back({n->pl().getX(), n->pl().getY()});
          break;
        }
      }
    }
  }

  double operator()(const GridNode* from, const std::set<GridNode*>& to) const {
    if (to.find(from->pl().getParent()) != to.end()) return 0;

    size_t ret = std::numeric_limits<size_t>::max();

    for (auto t : hull) {
      size_t temp = g->heurCost(from->pl().getX(), from->pl().getY(), t.first, t.second);
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
  TransitGraph draw(TransitGraph* g, GridGraph** gg, const Penalties& pens,
                    double gridSize, double borderRad);

 private:
  double getMaxDis(CombNode* to, CombEdge* e, double gridSize);
  void removeEdgesShorterThan(TransitGraph* g, double d);
  PolyLine<double> buildPolylineFromRes(const std::vector<GridEdge*>& l);
  double getCostFromRes(const std::vector<GridEdge*>& l);
  void addResidentEdges(gridgraph::GridGraph* g, CombEdge* e,
                        const std::vector<GridEdge*>& res);

  NodeCost writeNdCosts(GridNode* n, CombNode* origNode, CombEdge* e,
                        GridGraph* g);
  void settleRes(GridNode* startGridNd, GridNode* toGridNd, GridGraph* gg,
                 CombNode* from, CombNode* to, const GrEdgList& res,
                 CombEdge* e, size_t gen);
};

}  // namespace octi

#endif  // OCTI_OCTILINEARIZER_H_
