// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_CONVEXHULLOCTIGRIDGRAPH_H_
#define OCTI_BASEGRAPH_CONVEXHULLOCTIGRIDGRAPH_H_

#include "octi/basegraph/OctiGridGraph.h"

using util::geo::DPolygon;

namespace octi {
namespace basegraph {

class ConvexHullOctiGridGraph : public OctiGridGraph {
 public:
  using GridGraph::neigh;
  ConvexHullOctiGridGraph(const DPolygon& hull, const util::geo::DBox& bbox, double cellSize, double spacer,
                const Penalties& pens)
      : OctiGridGraph(bbox, cellSize, spacer, pens), _hull(hull) {
  }
  virtual void init();

 protected:
  virtual bool skip(size_t x, size_t y) const;
  virtual GridNode* writeNd(size_t x, size_t y);
  virtual GridNode* getNode(size_t x, size_t y) const;
 private:
  DPolygon _hull;
  std::vector<size_t> _ndIdx;
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_CONVEXHULLOCTIGRIDGRAPH_H_
