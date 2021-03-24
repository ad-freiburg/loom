// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_HEXGRIDGRAPH_H_
#define OCTI_BASEGRAPH_HEXGRIDGRAPH_H_

#include "octi/basegraph/GridGraph.h"

#define A 0.86602540378443865

namespace octi {
namespace basegraph {

class HexGridGraph : public GridGraph {
 public:
  using GridGraph::neigh;
  HexGridGraph(const util::geo::DBox& bbox, double cellSize, double spacer,
               const Penalties& pens)
      : GridGraph(bbox, cellSize, spacer, pens) {
    _a = _cellSize;
    _h = _a * A;
    _grid = Grid<GridNode*, Point, double>(_a, _h, bbox, false);

    _bendCosts[0] = _c.p_45 - _c.p_135;
    _bendCosts[2] = _c.p_45;
    _bendCosts[1] = _bendCosts[0] + _bendCosts[2];
  }

  virtual void init();
  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
  getHeur(const std::set<GridNode*>& to) const;
  virtual size_t maxDeg() const;
  virtual std::vector<double> getCosts() const;

 protected:
  virtual void writeInitialCosts();
  virtual double getBendPen(size_t i, size_t j) const;
  virtual GridNode* writeNd(size_t x, size_t y);
  virtual GridNode* neigh(size_t cx, size_t cy, size_t i) const;
  virtual GridNode* getNode(size_t x, size_t y) const;
  virtual size_t ang(size_t i, size_t j) const;

 private:
  double _a, _h;
  double _bendCosts[4];
};

struct HexGridGraphHeur
    : public util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float> {
  HexGridGraphHeur(const basegraph::GridGraph* g, const std::set<GridNode*>& to)
      : g(g), to(0) {}

  float operator()(const GridNode* from, const std::set<GridNode*>& to) const {
    return 0;
  }

  const octi::basegraph::BaseGraph* g;
  GridNode* to;
  std::vector<size_t> hull;
  float cheapestSink;
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_HEXGRIDGRAPH_H_
