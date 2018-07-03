// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_NODE_H_
#define UTIL_GRAPH_NODE_H_

#include <vector>

namespace util {
namespace graph {

// forward declaration of Edge
template <typename N, typename E>
class Edge;

template <typename N, typename E>
class Node {
 public:
  virtual const std::vector<Edge<N, E>*>& getAdjList() const = 0;
  virtual const std::vector<Edge<N, E>*>& getAdjListOut() const = 0;
  virtual const std::vector<Edge<N, E>*>& getAdjListIn() const = 0;

  virtual size_t getDeg() const = 0;
  virtual size_t getInDeg() const = 0;
  virtual size_t getOutDeg() const = 0;

  virtual bool hasEdgeIn(const Edge<N, E>* e) const = 0;
  virtual bool hasEdgeOut(const Edge<N, E>* e) const = 0;
  virtual bool hasEdge(const Edge<N, E>* e) const = 0;

  // add edge to this node's adjacency lists
  virtual void addEdge(Edge<N, E>* e) = 0;
  virtual void removeEdge(Edge<N, E>* e) = 0;

  virtual ~Node() {};

  virtual N& pl() = 0;
  virtual const N& pl() const = 0;
};

}}

#endif  // UTIL_GRAPH_NODE_H_
