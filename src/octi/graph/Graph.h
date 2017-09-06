// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/Graph.h"
#include "octi/graph/NodePL.h"
#include "octi/graph/EdgePL.h"

using util::graph::Graph;
using util::graph::Node;
using util::geo::Grid;
using util::geo::Point;

namespace octi {
namespace graph {

class Graph : public Graph<NodePL, EdgePL> {
 public:
  Graph();

 private:
};

}
}
