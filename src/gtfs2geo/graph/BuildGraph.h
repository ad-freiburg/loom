// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2GEO_GRAPH_BUILDGRAPH_H_
#define GTFS2GEO_GRAPH_BUILDGRAPH_H_


#include "util/graph/UndirGraph.h"

namespace gtfs2geo {
namespace graph {

class NodePL;
class EdgePL;

typedef util::graph::UndirGraph<NodePL, EdgePL> BuildGraph;
typedef util::graph::Node<NodePL, EdgePL> Node;
typedef util::graph::Edge<NodePL, EdgePL> Edge;

}
}

#endif  // GTFS2GEO_GRAPH_BUILDGRAPH_H_
