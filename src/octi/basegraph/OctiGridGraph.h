// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_OCTIGRIDGRAPH_H_
#define OCTI_BASEGRAPH_OCTIGRIDGRAPH_H_

#include "octi/basegraph/GridGraph.h"

namespace octi {
namespace basegraph {

class OctiGridGraph : public GridGraph {
 public:
  OctiGridGraph(const util::geo::DBox& bbox, double cellSize, double spacer,
                const Penalties& pens)
      : GridGraph(bbox, cellSize, spacer, pens) {}
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_OCTIGRIDGRAPH_H_
