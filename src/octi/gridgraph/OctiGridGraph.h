// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_OCTIGRIDGRAPH_H_
#define OCTI_GRIDGRAPH_OCTIGRIDGRAPH_H_

#include "octi/gridgraph/GridGraph.h"

namespace octi {
namespace gridgraph {

class OctiGridGraph : public GridGraph {
 public:
  OctiGridGraph(const util::geo::DBox& bbox, double cellSize, double spacer,
                const Penalties& pens)
      : GridGraph(bbox, cellSize, spacer, pens) {}
};
}  // namespace gridgraph
}  // namespace octi

#endif  // OCTI_GRIDGRAPH_OCTIGRIDGRAPH_H_
