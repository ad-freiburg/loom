// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_GRAPH_H_
#define UTIL_GRAPH_GRAPH_H_

#include <set>
#include <string>

#include "util/graph/Edge.h"
#include "util/graph/Node.h"

namespace util {
namespace graph {

template <typename N, typename E>
class Graph {
 public:
  explicit Graph();
  explicit Graph(bool directed);
  ~Graph();

  Node<N, E>* addNode(Node<N, E>* n);
  Edge<N, E>* addEdge(Node<N, E>* from, Node<N, E>* to, const E& p);
  Edge<N, E>* getEdge(Node<N, E>* from, Node<N, E>* to);

  void deleteNode(Node<N, E>* n);
  void deleteEdge(Node<N, E>* from, Node<N, E>* to);

  const std::set<Node<N, E>*>& getNodes() const;
  std::set<Node<N, E>*>* getNodes();

  Node<N, E>* mergeNodes(Node<N, E>* a, Node<N, E>* b);

 private:
  std::set<Node<N, E>*> _nodes;
  bool _directed;
};

#include "util/graph/Graph.tpp"
}
}

#endif  // UTIL_GRAPH_GRAPH_H_
