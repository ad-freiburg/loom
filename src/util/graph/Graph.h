// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_GRAPH_H_
#define UTIL_GRAPH_GRAPH_H_

#include <set>
#include <string>
#include <iostream>
#include <cassert>

#include "util/graph/Edge.h"
#include "util/graph/Node.h"

namespace util {
namespace graph {

template <typename N, typename E>
class Graph {
 public:
  virtual ~Graph();
  virtual Node<N, E>* addNd() = 0;
  virtual Node<N, E>* addNd(const N& pl) = 0;
  Edge<N, E>* addEdg(Node<N, E>* from, Node<N, E>* to);
  virtual Edge<N, E>* addEdg(Node<N, E>* from, Node<N, E>* to, const E& p) = 0;
  Edge<N, E>* getEdg(Node<N, E>* from, Node<N, E>* to);
  const Edge<N, E>* getEdg(Node<N, E>* from, Node<N, E>* to) const;

  virtual Node<N, E>* mergeNds(Node<N, E>* a, Node<N, E>* b) = 0;

  const std::set<Node<N, E>*>& getNds() const;
  std::set<Node<N, E>*>* getNds();

  typename std::set<Node<N, E>*>::iterator delNd(Node<N, E>* n);
  typename std::set<Node<N, E>*>::iterator delNd(
      typename std::set<Node<N, E>*>::iterator i);
  void delEdg(Node<N, E>* from, Node<N, E>* to);

 private:
  std::set<Node<N, E>*> _nodes;
};

#include "util/graph/Graph.tpp"
}
}

#endif  // UTIL_GRAPH_GRAPH_H_
