// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_OCTIHANANGRAPH_H_
#define OCTI_BASEGRAPH_OCTIHANANGRAPH_H_

#include "octi/basegraph/OctiGridGraph.h"

namespace octi {
namespace basegraph {

class OctiHananGraph : public OctiGridGraph {
 public:
  using OctiGridGraph::neigh;
  OctiHananGraph(const util::geo::DBox& bbox, const combgraph::CombGraph& cg,
                 double cellSize, double spacer, size_t iters,
                 const Penalties& pens)
      : OctiGridGraph(bbox, cellSize, spacer, pens), _cg(cg), _iters(iters) {}

  virtual void unSettleEdg(CombEdge* ce, GridNode* a, GridNode* b);
  virtual void settleEdg(GridNode* a, GridNode* b, CombEdge* e);
  virtual CrossEdgPairs getCrossEdgPairs() const;
  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual size_t maxDeg() const;
  virtual double ndMovePen(const CombNode* cbNd, const GridNode* grNd) const;
  virtual void init();

 protected:
  virtual GridNode* writeNd(size_t x, size_t y);
  virtual GridNode* neigh(size_t cx, size_t cy, size_t i) const;
  virtual GridNode* getNode(size_t x, size_t y) const;
  virtual double getBendPen(size_t i, size_t j) const;
  virtual size_t ang(size_t i, size_t j) const;
  virtual void connectNodes(GridNode* grNdA, GridNode* grNdB, size_t dir);
  virtual void writeInitialCosts();
  std::set<std::pair<size_t, size_t>> getIterCoords(
      const std::set<std::pair<size_t, size_t>>& inCoords) const;

  const combgraph::CombGraph& _cg;
  size_t _iters;
  std::vector<size_t> _ndIdx;
  std::vector<GridNode*> _neighs;

  std::map<GridEdge*, std::vector<std::pair<GridEdge*, GridEdge*>>> _edgePairs;
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_OCTIHANANGRAPH_H_
