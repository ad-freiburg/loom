// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_EDGE_H_
#define UTIL_GRAPH_EDGE_H_

#include <vector>
#include "util/graph/Node.h"

namespace util {
namespace graph {

template <typename N, typename E>
class Edge {
 public:
  Edge(Node<N, E>* from, Node<N, E>* to, const E& pl);

  Node<N, E>* getFrom() const;
  Node<N, E>* getTo() const;

  Node<N, E>* getOtherNd(const Node<N, E>* notNode) const;

  E& pl();
  const E& pl() const;

 private:
  Node<N, E>* _from;
  Node<N, E>* _to;
  E _pl;
};

#include "util/graph/Edge.tpp"

}}

#endif  // UTIL_GRAPH_EDGE_H_

