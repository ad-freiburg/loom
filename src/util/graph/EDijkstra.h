// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_EDIJKSTRA_H_
#define UTIL_GRAPH_EDIJKSTRA_H_

#include <limits>
#include <list>
#include <queue>
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

// edge-based dijkstra - settles edges instead of nodes
class EDijkstra : public ShortestPath<EDijkstra> {
 public:
  template <typename N, typename E, typename C>
  struct RouteEdge {
    RouteEdge() : e(0), parent(0), d(), h(), n(0) {}
    RouteEdge(Edge<N, E>* e) : e(e), parent(0), d(), h(), n(0) {}
    RouteEdge(Edge<N, E>* e, Edge<N, E>* parent, Node<N, E>* n, C d)
        : e(e), parent(parent), d(d), h(), n(n) {}
    RouteEdge(Edge<N, E>* e, Edge<N, E>* parent, Node<N, E>* n, C d, C h)
        : e(e), parent(parent), d(d), h(h), n(n) {}

    Edge<N, E>* e;
    Edge<N, E>* parent;

    C d;
    C h;

    Node<N, E>* n;

    bool operator<(const RouteEdge<N, E, C>& p) const {
      return h > p.h || (h == p.h && d > p.d);
    }
  };

  template <typename N, typename E, typename C>
  struct CostFunc : public ShortestPath::CostFunc<N, E, C> {
    C operator()(const Node<N, E>* from, const Edge<N, E>* e,
                 const Node<N, E>* to) const {
      UNUSED(from);
      UNUSED(e);
      UNUSED(to);
      return C();
    };
  };

  template <typename N, typename E, typename C>
  struct HeurFunc : public ShortestPath::HeurFunc<N, E, C> {
    C operator()(const Node<N, E>* from,
                 const std::set<Node<N, E>*>& to) const {
      UNUSED(from);
      UNUSED(to);
      return C();
    };
  };

  template <typename N, typename E, typename C>
  using Settled = std::unordered_map<Edge<N, E>*, RouteEdge<N, E, C> >;

  template <typename N, typename E, typename C>
  using PQ = std::priority_queue<RouteEdge<N, E, C> >;

  template <typename N, typename E, typename C>
  static C shortestPathImpl(const std::set<Edge<N, E>*> from,
                            const std::set<Edge<N, E>*>& to,
                            const ShortestPath::CostFunc<N, E, C>& costFunc,
                            const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                            const ShortestPath::CostFunc<N, E, C>& costFunc,
                            const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(Edge<N, E>* from, const std::set<Node<N, E>*>& to,
                            const ShortestPath::CostFunc<N, E, C>& costFunc,
                            const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(const std::set<Edge<N, E>*>& from,
                            const std::set<Node<N, E>*>& to,
                            const ShortestPath::CostFunc<N, E, C>& costFunc,
                            const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathImpl(
      const std::set<Edge<N, E>*>& from,
      const ShortestPath::CostFunc<N, E, C>& costFunc, bool rev);

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathImpl(
      Edge<N, E>* from, const std::set<Edge<N, E>*>& to,
      const ShortestPath::CostFunc<N, E, C>& costFunc,
      const ShortestPath::HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes);

  template <typename N, typename E, typename C>
  static void buildPath(Edge<N, E>* curE, const Settled<N, E, C>& settled,
                        NList<N, E>* resNodes, EList<N, E>* resEdges);

  template <typename N, typename E, typename C>
  static inline void relax(RouteEdge<N, E, C>& cur,
                           const std::set<Edge<N, E>*>& to,
                           const ShortestPath::CostFunc<N, E, C>& costFunc,
                           const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                           PQ<N, E, C>& pq);

  template <typename N, typename E, typename C>
  static void relaxInv(RouteEdge<N, E, C>& cur,
                      const ShortestPath::CostFunc<N, E, C>& costFunc,
                      PQ<N, E, C>& pq);

  static size_t ITERS;
};

#include "util/graph/EDijkstra.tpp"
}
}

#endif  // UTIL_GRAPH_DIJKSTRA_H_
