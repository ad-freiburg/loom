// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/Geo.h"
#include "util/graph/Graph.h"
#include "octi/gridgraph/NodePL.h"
#include "octi/gridgraph/EdgePL.h"

using util::graph::Graph;

namespace octi {
namespace gridgraph {

class GridGraph : public Graph<NodePL, EdgePL> {
 public:
  GridGraph(util::geo::Box bbox, double cellSize);

 private:
  util::geo::Box _bbox;
  double _cellSize;
};

}
}
