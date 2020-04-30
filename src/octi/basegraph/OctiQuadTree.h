// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_OCTIQUADTREE_H_
#define OCTI_BASEGRAPH_OCTIQUADTREE_H_

#include "octi/basegraph/OctiGridGraph.h"

namespace octi {
namespace basegraph {

class OctiQuadTree : public OctiGridGraph {
 public:
  using OctiGridGraph::neigh;
  OctiQuadTree(const util::geo::DBox& bbox, const combgraph::CombGraph& cg,
               double cellSize, double spacer, const Penalties& pens)
      : OctiGridGraph(bbox, cellSize, spacer, pens), _cg(cg) {}

  virtual void unSettleEdg(GridNode* a, GridNode* b);
  virtual void settleEdg(GridNode* a, GridNode* b, CombEdge* e, size_t order);
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
  virtual double heurCost(int64_t xa, int64_t ya, int64_t xb, int64_t yb) const;
  virtual size_t getGrNdDeg(const CombNode* nd, size_t x, size_t y) const;

 private:
  void connectNodes(GridNode* grNdA, GridNode* grNdB, size_t dir);

  const combgraph::CombGraph& _cg;
  std::map<std::pair<size_t, size_t>, size_t> _ndIdx;
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_OCTIQUADTREE_H_
