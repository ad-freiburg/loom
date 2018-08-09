// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_DIJKSTRA_H_
#define UTIL_GRAPH_DIJKSTRA_H_

#include <limits>
#include <list>
#include <queue>
#include <set>
#include <set>
#include <unordered_map>
#include "util/graph/Edge.h"
#include "util/graph/Graph.h"
#include "util/graph/Node.h"
#include "util/graph/ShortestPath.h"

namespace util {
namespace graph {

using util::graph::Graph;
using util::graph::Node;
using util::graph::Edge;

// dijkstras algorithm for util graph
class Dijkstra : public ShortestPath<Dijkstra> {
 public:
  template <typename N, typename E, typename C>
  struct RouteNode {
    RouteNode() : n(0), parent(0), d(), h(), e(0) {}
    RouteNode(Node<N, E>* n) : n(n), parent(0), d(), h(), e(0) {}
    RouteNode(Node<N, E>* n, Node<N, E>* parent, C d, Edge<N, E>* e)
        : n(n), parent(parent), d(d), h(), e(e) {}
    RouteNode(Node<N, E>* n, Node<N, E>* parent, C d, C h, Edge<N, E>* e)
        : n(n), parent(parent), d(d), h(h), e(e) {}

    Node<N, E>* n;
    Node<N, E>* parent;

    C d;
    C h;

    Edge<N, E>* e;

    bool operator<(const RouteNode<N, E, C>& p) const {
      return h > p.h || (h == p.h && d > p.d);
    }
  };

  template <typename N, typename E, typename C>
  using Settled = std::unordered_map<Node<N, E>*, RouteNode<N, E, C> >;

  template <typename N, typename E, typename C>
  using PQ = std::priority_queue<RouteNode<N, E, C> >;

  template <typename N, typename E, typename C>
  struct CostFunc : public ShortestPath::CostFunc<N, E, C> {
    C operator()(const Edge<N, E>* from, const Node<N, E>* n,
                 const Edge<N, E>* to) const {
      UNUSED(from);
      UNUSED(n);
      UNUSED(to);
      return C();
    };
  };

  template <typename N, typename E, typename C>
  struct HeurFunc : public ShortestPath::HeurFunc<N, E, C> {
    C operator()(const Edge<N, E>* from,
                 const std::set<Edge<N, E>*>& to) const {
      UNUSED(from);
      UNUSED(to);
      return C();
    };
  };

  template <typename N, typename E, typename C>
  static std::unordered_map<Node<N, E>*, C> shortestPathImpl(
      Node<N, E>* from, const std::set<Node<N, E>*>& to,
      const ShortestPath::CostFunc<N, E, C>& costFunc, const ShortestPath::HeurFunc<N, E, C>&,
      std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Node<N, E>*, NList<N, E>*> resNode);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(const std::set<Node<N, E>*> from, const std::set<Node<N, E>*>& to,
                            const ShortestPath::CostFunc<N, E, C>& costFunc,
                            const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                            const ShortestPath::CostFunc<N, E, C>& costFunc,
                            const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static void relax(RouteNode<N, E, C>& cur, const std::set<Node<N, E>*>& to,
                    const ShortestPath::CostFunc<N, E, C>& costFunc,
                    const ShortestPath::HeurFunc<N, E, C>& heurFunc, PQ<N, E, C>& pq);

  template <typename N, typename E, typename C>
  static void buildPath(Node<N, E>* curN, Settled<N, E, C>& settled,
                        NList<N, E>* resNodes, EList<N, E>* resEdges);

  static size_t ITERS;
};

#include "util/graph/Dijkstra.tpp"
}
}

#endif  // UTIL_GRAPH_DIJKSTRA_H_
