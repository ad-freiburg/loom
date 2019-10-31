// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_OCTILINEARIZER_H_
#define OCTI_OCTILINEARIZER_H_

#include <unordered_set>
#include <vector>
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
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

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdge;

using octi::combgraph::EdgeOrdering;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::CombEdge;
using octi::combgraph::Drawing;

using util::graph::Dijkstra;

typedef Dijkstra::EList<GridNodePL, GridEdgePL> GrEdgList;
typedef Dijkstra::NList<GridNodePL, GridEdgePL> GrNdList;
typedef std::pair<std::set<GridNode*>, std::set<GridNode*>> RtPair;
typedef std::map<CombNode*, std::pair<size_t, size_t>> SettledPos;

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
  GridHeur(GridGraph* g, const std::set<GridNode*>& to) : g(g), to(0) {
    if (to.size() == 1) this->to = *to.begin();

    cheapestSink = std::numeric_limits<double>::infinity();

    for (auto n : to) {
      size_t i = 0;
      for (; i < 8; i++) {
        auto neigh = g->getNeighbor(n->pl().getX(), n->pl().getY(), i);
        if (neigh && to.find(neigh) == to.end()) {
          hull.push_back({n->pl().getX(), n->pl().getY()});
          break;
        }
      }
      for (size_t j = i; j < 8; j++) {
        double sinkCost = g->getEdg(n, n->pl().getPort(j))->pl().cost();
        if (sinkCost < cheapestSink) cheapestSink = sinkCost;
      }
    }
  }

  double operator()(const GridNode* from, const std::set<GridNode*>& to) const {
    if (to.find(from->pl().getParent()) != to.end()) return 0;

    double ret = std::numeric_limits<double>::infinity();

    for (auto t : hull) {
      double temp =
          g->heurCost(from->pl().getParent()->pl().getX(),
                      from->pl().getParent()->pl().getY(), t.first, t.second);
      if (temp < ret) ret = temp;
    }

    return ret + cheapestSink;
  }

  octi::gridgraph::GridGraph* g;
  GridNode* to;
  std::vector<std::pair<size_t, size_t>> hull;
  double cheapestSink;
};

class Octilinearizer {
 public:
  Octilinearizer() {}
  TransitGraph draw(TransitGraph* g, GridGraph** gg, const Penalties& pens,
                    double gridSize, double borderRad);

  TransitGraph drawILP(TransitGraph* g, GridGraph** gg, const Penalties& pens,
                       double gridSize, double borderRad);

 private:
  void removeEdgesShorterThan(TransitGraph* g, double d);

  void writeNdCosts(GridNode* n, CombNode* origNode, CombEdge* e, GridGraph* g);

  void settleRes(GridNode* startGridNd, GridNode* toGridNd, GridGraph* gg,
                 CombNode* from, CombNode* to, const GrEdgList& res,
                 CombEdge* e);

  std::vector<CombEdge*> getOrdering(const CombGraph& cg) const;

  bool draw(const std::vector<CombEdge*>& order, GridGraph* gg,
            Drawing* drawing);
  bool draw(const std::vector<CombEdge*>& order, const SettledPos& settled,
            GridGraph* gg, Drawing* drawing);

  SettledPos getNeighbor(const SettledPos& pos, const std::vector<CombNode*>&,
                         size_t i) const;

  RtPair getRtPair(CombNode* frCmbNd, CombNode* toCmbNd,
                   const SettledPos& settled, GridGraph* gg);

  std::set<GridNode*> getCands(CombNode* cmBnd, const SettledPos& settled,
                               GridGraph* gg, double maxDis);
};

}  // namespace octi

#endif  // OCTI_OCTILINEARIZER_H_
