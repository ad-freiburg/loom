// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_ALGORITHM_H_
#define UTIL_GRAPH_ALGORITHM_H_

#include <stack>
#include "util/graph/Edge.h"
#include "util/graph/Node.h"
#include "util/graph/UndirGraph.h"

namespace util {
namespace graph {

using util::graph::Graph;
using util::graph::Node;
using util::graph::Edge;

// collection of general graph algorithms
class Algorithm {
 public:
  template <typename N, typename E>
  struct EdgeCheckFunc {
    virtual bool operator()(const Node<N, E>* frNd,
                            const Edge<N, E>* edge) const {
      UNUSED(frNd);
      UNUSED(edge);
      return true;
    };
  };

  template <typename N, typename E>
  static std::vector<std::set<Node<N, E>*> > connectedComponents(
      const UndirGraph<N, E>& g);

  template <typename N, typename E>
  static std::vector<std::set<Node<N, E>*> > connectedComponents(
      const UndirGraph<N, E>& g, const EdgeCheckFunc<N, E>& checkFunc);
};

#include "util/graph/Algorithm.tpp"
}
}

#endif  // UTIL_GRAPH_ALGORITHM_H_
