// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_GRIDGRAPH_H_
#define OCTI_BASEGRAPH_GRIDGRAPH_H_

#include <queue>
#include <set>
#include <unordered_map>
#include "octi/basegraph/BaseGraph.h"
#include "octi/basegraph/GridEdgePL.h"
#include "octi/basegraph/GridNodePL.h"
#include "octi/basegraph/NodeCost.h"
#include "octi/combgraph/CombGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/Dijkstra.h"
#include "util/graph/DirGraph.h"

#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"

using util::geo::Grid;
using util::geo::Point;
using util::graph::DirGraph;
using util::graph::Node;

using octi::combgraph::CombEdge;
using octi::combgraph::CombNode;

namespace octi {
namespace basegraph {

class GridGraph : public BaseGraph {
 public:
  GridGraph(const util::geo::DBox& bbox, double cellSize, double spacer,
            const Penalties& pens);

  virtual double getCellSize() const;

  virtual NodeCost nodeBendPen(GridNode* n, CombNode* origNd, CombEdge* e);
  virtual NodeCost topoBlockPen(GridNode* n, CombNode* origNode, CombEdge* e);
  virtual NodeCost spacingPen(GridNode* n, CombNode* origNode, CombEdge* e);

  virtual double heurCost(int64_t xa, int64_t ya, int64_t xb, int64_t yb) const;

  virtual std::priority_queue<Candidate> getGridNdCands(
      const util::geo::DPoint& p, size_t maxGrD) const;

  virtual void addCostVec(GridNode* n, const NodeCost& addC);

  virtual void openSinkTo(GridNode* n, double cost);
  virtual void closeSinkTo(GridNode* n);
  virtual void openSinkFr(GridNode* n, double cost);
  virtual void closeSinkFr(GridNode* n);
  virtual void openTurns(GridNode* n);
  virtual void closeTurns(GridNode* n);

  virtual GridNode* neigh(const GridNode* n, size_t i) const;
  virtual size_t maxDeg() const;

  virtual std::set<GridNode*> getGrNdCands(CombNode* n, size_t maxGrDist);

  virtual void settleNd(GridNode* n, CombNode* cn);
  virtual void settleEdg(GridNode* a, GridNode* b, CombEdge* e, size_t rndrOr);

  virtual const Penalties& getPens() const;

  virtual void unSettleEdg(CombEdge* ce, GridNode* a, GridNode* b);
  virtual void unSettleNd(CombNode* a);

  virtual bool isSettled(const CombNode* cn);

  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual void init();
  virtual void reset();

  virtual GridNode* getSettled(const CombNode* cnd) const;

  virtual double ndMovePen(const CombNode* cbNd, const GridNode* grNd) const;

  virtual GridNode* getGrNdById(size_t id) const;
  virtual const GridEdge* getGrEdgById(std::pair<size_t, size_t> id) const;
  virtual void addResEdg(GridEdge* ge, CombEdge* cg);
  virtual std::set<CombEdge*> getResEdgs(GridEdge* ge);

  virtual CrossEdgPairs getCrossEdgPairs() const;

  virtual void writeGeoCoursePens(const CombEdge* ce, GeoPensMap* target,
                                  double pen);

  virtual void addObstacle(const util::geo::Polygon<double>& obst);

  virtual const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
  getHeur(const std::set<GridNode*>& to) const;

  virtual PolyLine<double> geomFromPath(
      const std::vector<std::pair<size_t, size_t>>& res) const;

  virtual void prunePorts();

 protected:
  util::geo::DBox _bbox;
  Penalties _c;

  Grid<GridNode*, Point, double> _grid;
  double _cellSize, _spacer;
  std::unordered_map<const CombNode*, GridNode*> _settled;

  double _heurHopCost;

  // encoding portable IDs for each node
  std::vector<GridNode*> _nds;

  // edge id counter
  size_t _edgeCount;

  std::vector<util::geo::Polygon<double>> _obstacles;

  // may be multiple resident edges if hard constraints are relaxed
  std::unordered_map<GridEdge*, std::set<CombEdge*>> _resEdgs;

  const Grid<GridNode*, Point, double>& getGrid() const;

  virtual void writeInitialCosts();
  virtual void writeObstacleCost(const util::geo::Polygon<double>& obst);
  virtual void reWriteObstCosts();

  virtual double getBendPen(size_t origI, size_t targetI) const;
  virtual size_t ang(size_t i, size_t j) const;

  virtual GridNode* getNode(size_t x, size_t y) const;

  virtual GridNode* writeNd(size_t x, size_t y);

  virtual GridNode* neigh(size_t cx, size_t cy, size_t i) const;

  virtual void getSettledAdjEdgs(GridNode* n, CombNode* origNd,
                                 CombEdge* outgoing[8]);

  virtual size_t getGrNdDeg(const CombNode* nd, size_t x, size_t y) const;

 private:
  double _bendCosts[2];
};

struct GridCost
    : public util::graph::Dijkstra::CostFunc<GridNodePL, GridEdgePL, float> {
  GridCost(float inf) : _inf(inf) {}
  virtual float operator()(const GridNode* from, const GridEdge* e,
                           const GridNode* to) const {
    UNUSED(from);
    UNUSED(to);
    // const_cast<GridNode*>(from)->pl().visited = true;
    return e->pl().cost();
  }

  float _inf;

  virtual float inf() const { return _inf; }
};

struct GridGraphHeur
    : public util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float> {
  GridGraphHeur(const basegraph::GridGraph* g, const std::set<GridNode*>& to)
      : g(g) {
    cheapestSink = std::numeric_limits<float>::infinity();

    for (auto n : to) {
      assert(n->pl().getParent() == n);
      size_t i = 0;
      for (; i < g->maxDeg(); i++) {
        if (!n->pl().getPort(i)) continue;
        float sinkCost = g->getEdg(n->pl().getPort(i), n)->pl().cost();
        if (sinkCost < cheapestSink) cheapestSink = sinkCost;
        auto neigh = g->neigh(n, i);
        if (neigh && to.find(neigh) == to.end()) {
          hull.push_back(n->pl().getX());
          hull.push_back(n->pl().getY());
          break;
        }
      }
      for (size_t j = i; j < g->maxDeg(); j++) {
        if (!n->pl().getPort(j)) continue;
        float sinkCost = g->getEdg(n->pl().getPort(j), n)->pl().cost();
        if (sinkCost < cheapestSink) cheapestSink = sinkCost;
      }
    }
  }

  float operator()(const GridNode* from, const std::set<GridNode*>& to) const {
    // const_cast<GridNode*>(from)->pl().visited = true;
    if (to.count(from->pl().getParent())) return 0;

    float ret = std::numeric_limits<float>::infinity();

    for (size_t i = 0; i < hull.size(); i += 2) {
      float tmp = g->heurCost(from->pl().getParent()->pl().getX(),
                              from->pl().getParent()->pl().getY(), hull[i],
                              hull[i + 1]);
      if (tmp < ret) ret = tmp;
    }

    return ret + cheapestSink;
  }

  const octi::basegraph::BaseGraph* g;
  std::vector<size_t> hull;
  float cheapestSink;
};

}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_GRIDGRAPH_H_
