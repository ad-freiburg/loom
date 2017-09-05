// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_DIJKSTRA_H_
#define UTIL_GRAPH_DIJKSTRA_H_

#include <unordered_map>
#include <queue>
#include <list>
#include "util/graph/Edge.h"
#include "util/graph/Graph.h"
#include "util/graph/Node.h"

using util::graph::Graph;
using util::graph::Node;
using util::graph::Edge;

namespace util {
namespace graph {

template <typename N, typename E>
struct RouteNode {
  RouteNode() : n(0), parent(0), d(0), e(0) {}
  RouteNode(Node<N, E>* n, Node<N, E>* parent, double d, Edge<N, E>* e)
      : n(n), parent(parent), d(d), e(e) {}

  Node<N, E>* n;
  Node<N, E>* parent;

  double d;
  Edge<N, E>* e;

  bool operator<(const RouteNode<N, E>& p) const { return d > p.d; }
};

class Dijkstra {
 public:
  template <typename N, typename E>
  static bool shortestPath(Node<N, E>* from, Node<N, E>* to,
                           std::list<Edge<N, E>*>* res);
};

#include "util/graph/Dijkstra.tpp"
}
}

#endif  // UTIL_GRAPH_DIJKSTRA_H_
