// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_PSEUDOORTHORADIALGRAPH_H_
#define OCTI_BASEGRAPH_PSEUDOORTHORADIALGRAPH_H_

#include "octi/basegraph/GridGraph.h"

namespace octi {
namespace basegraph {

class PseudoOrthoRadialGraph : public GridGraph {
 public:
  using GridGraph::neigh;
  PseudoOrthoRadialGraph(const util::geo::DBox& bbox, double cellSize,
                         double spacer, const Penalties& pens)
      : GridGraph(util::geo::extendBox(
                      bbox, util::geo::getBoundingBox(util::geo::rotate(
                                convexHull(bbox), 90, centroid(bbox)))),
                  cellSize, spacer, pens) {
    double circum = cellSize * 1.3 * 2 * M_PI;
    _numBeams = circum / cellSize;
  }

  virtual void init();
  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
  getHeur(const std::set<GridNode*>& to) const;

  virtual PolyLine<double> geomFromPath(
      const std::vector<std::pair<size_t, size_t>>& res) const;
  virtual double ndMovePen(const CombNode* cbNd, const GridNode* grNd) const;

 protected:
  virtual void writeInitialCosts();
  virtual GridNode* writeNd(size_t x, size_t y);
  virtual GridNode* neigh(size_t cx, size_t cy, size_t i) const;
  virtual GridNode* getNode(size_t x, size_t y) const;

 private:
  virtual int multi(size_t y) const;

  size_t _numBeams;
};

struct PseudoOrthoRadialGraphHeur
    : public util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float> {
  PseudoOrthoRadialGraphHeur(const basegraph::GridGraph* g,
                             const std::set<GridNode*>& to)
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

#endif  // OCTI_BASEGRAPH_PSEUDOORTHORADIALGRAPH_H_
