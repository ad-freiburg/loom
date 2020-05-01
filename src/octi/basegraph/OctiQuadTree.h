// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_OCTIQUADTREE_H_
#define OCTI_BASEGRAPH_OCTIQUADTREE_H_

#include "octi/basegraph/OctiHananGraph.h"

namespace octi {
namespace basegraph {

class OctiQuadTree : public OctiHananGraph {
 public:
  using OctiGridGraph::neigh;
  OctiQuadTree(const util::geo::DBox& bbox, const combgraph::CombGraph& cg,
               double cellSize, double spacer, const Penalties& pens)
      : OctiHananGraph(bbox, cg, cellSize, spacer, pens) {}

  virtual void unSettleEdg(GridNode* a, GridNode* b);
  virtual void settleEdg(GridNode* a, GridNode* b, CombEdge* e, size_t order);
  virtual CrossEdgPairs getCrossEdgPairs() const;
  virtual double ndMovePen(const CombNode* cbNd, const GridNode* grNd) const;
  virtual void init();
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_OCTIQUADTREE_H_
