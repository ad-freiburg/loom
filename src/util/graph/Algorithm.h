// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_ALGORITHM_H_
#define UTIL_GRAPH_ALGORITHM_H_

#include <stack>
#include "util/graph/Edge.h"
#include "util/graph/UndirGraph.h"
#include "util/graph/Node.h"

namespace util {
namespace graph {

using util::graph::Graph;
using util::graph::Node;
using util::graph::Edge;

// collection of general graph algorithms
class Algorithm {
 public:
  template <typename N, typename E>
  static std::vector<std::set<Node<N, E>*> > connectedComponents(
      const UndirGraph<N, E>& g);
};

#include "util/graph/Algorithm.tpp"

}
}

#endif  // UTIL_GRAPH_ALGORITHM_H_
