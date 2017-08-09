// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_NODE_H_
#define UTIL_GRAPH_NODE_H_

#include <set>
#include "transitmap/graph/Penalties.h"

namespace util {
namespace graph {

// forward declaration of Edge
template <typename N, typename E>
class Edge;

template <typename N, typename E>
class Node {
 public:
  Node(const N& pl);
  ~Node();

  const std::set<Edge<N, E>*>& getAdjListOut() const;
  const std::set<Edge<N, E>*>& getAdjListIn() const;

  std::set<Edge<N, E>*> getAdjList() const;

  // add edge to this node's adjacency lists
  void addEdge(Edge<N, E>* e);

  // remove edge from this node's adjacency lists
  void removeEdge(Edge<N, E>* e);

  N& pl();
  const N& pl() const;

 private:
  std::set<Edge<N, E>*> _adjListIn;
  std::set<Edge<N, E>*> _adjListOut;
  N _pl;
};

#include "util/graph/Node.tpp"

}}

#endif  // UTIL_GRAPH_NODE_H_
