// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_ORTHORADIALGRAPH_H_
#define OCTI_BASEGRAPH_ORTHORADIALGRAPH_H_

#include "octi/basegraph/GridGraph.h"

namespace octi {
namespace basegraph {

class OrthoRadialGraph : public GridGraph {
 public:
  using GridGraph::getNeighbor;
  OrthoRadialGraph(size_t numBeams, const util::geo::DBox& bbox,
                   double cellSize, double spacer, const Penalties& pens)
      : GridGraph(bbox, cellSize, spacer, pens), _numBeams(numBeams) {}

  virtual void init();
  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
  getHeur(const std::set<GridNode*>& to) const;

 protected:
  virtual void writeInitialCosts();
  virtual GridNode* writeNd(size_t x, size_t y);
  virtual GridNode* getNeighbor(size_t cx, size_t cy, size_t i) const;
  virtual GridNode* getNode(size_t x, size_t y) const;

 private:
  size_t _numBeams;
};

struct OrthoRadialGraphHeur
    : public util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float> {
  OrthoRadialGraphHeur(const basegraph::GridGraph* g,
                       const std::set<GridNode*>& to)
      : g(g), to(0) {
  }

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

#endif  // OCTI_BASEGRAPH_ORTHORADIALGRAPH_H_
