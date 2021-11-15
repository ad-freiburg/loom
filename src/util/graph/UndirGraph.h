// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_UNDIRGRAPH_H_
#define UTIL_GRAPH_UNDIRGRAPH_H_

#include <set>
#include <string>

#include "util/graph/Graph.h"
#include "util/graph/Edge.h"
#include "util/graph/UndirNode.h"

namespace util {
namespace graph {

template <typename N, typename E>
using UndirEdge = Edge<N, E>;

template <typename N, typename E>
class UndirGraph : public Graph<N, E> {
 public:
  explicit UndirGraph();

  using Graph<N, E>::addEdg;

  Node<N, E>* addNd();
  Node<N, E>* addNd(UndirNode<N, E>* n);
  Node<N, E>* addNd(const N& pl);
  Edge<N, E>* addEdg(Node<N, E>* from, Node<N, E>* to, const E& p);

  virtual Node<N, E>* mergeNds(Node<N, E>* a, Node<N, E>* b);

};

#include "util/graph/UndirGraph.tpp"
}
}

#endif  // UTIL_GRAPH_UNDIRGRAPH_H_
