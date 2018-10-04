// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_DIRNODE_H_
#define UTIL_GRAPH_DIRNODE_H_

#include <vector>
#include <algorithm>
#include "util/graph/Node.h"

namespace util {
namespace graph {

// forward declaration of Edge

template <typename N, typename E>
class DirNode : public Node<N, E> {
 public:
  DirNode();
  DirNode(const N& pl);
  ~DirNode();

  const std::vector<Edge<N, E>*>& getAdjList() const;
  const std::vector<Edge<N, E>*>& getAdjListIn() const;
  const std::vector<Edge<N, E>*>& getAdjListOut() const;

  size_t getDeg() const;
  size_t getInDeg() const;
  size_t getOutDeg() const;

  bool hasEdgeIn(const Edge<N, E>* e) const;
  bool hasEdgeOut(const Edge<N, E>* e) const;
  bool hasEdge(const Edge<N, E>* e) const;

  // add edge to this node's adjacency lists
  void addEdge(Edge<N, E>* e);

  // remove edge from this node's adjacency lists
  void removeEdge(Edge<N, E>* e);

  N& pl();
  const N& pl() const;

 private:
  std::vector<Edge<N, E>*> _adjListIn;
  std::vector<Edge<N, E>*> _adjListOut;
  N _pl;

  bool adjInContains(const Edge<N, E>* e) const;
  bool adjOutContains(const Edge<N, E>* e) const;
};

#include "util/graph/DirNode.tpp"

}}

#endif  // UTIL_GRAPH_DIRNODE_H_
