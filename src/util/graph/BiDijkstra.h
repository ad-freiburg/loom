// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_BIDIJKSTRA_H_
#define UTIL_GRAPH_BIDIJKSTRA_H_

#include <limits>
#include <list>
#include <queue>
#include <set>
#include <algorithm>
#include <unordered_map>
#include "util/graph/Edge.h"
#include "util/graph/Graph.h"
#include "util/graph/Node.h"
#include "util/graph/ShortestPath.h"

namespace util {
namespace graph {

using util::graph::Edge;
using util::graph::Graph;
using util::graph::Node;

// bidirectional dijkstras algorithm for util graph
class BiDijkstra : public ShortestPath<BiDijkstra> {
 public:
  template <typename N, typename E, typename C>
  struct RouteNode {
    RouteNode() : n(0), parent(0), d(), h() {}
    RouteNode(Node<N, E>* n) : n(n), parent(0), d(), h() {}
    RouteNode(Node<N, E>* n, Node<N, E>* parent, C d)
        : n(n), parent(parent), d(d), h() {}
    RouteNode(Node<N, E>* n, Node<N, E>* parent, C d, C h)
        : n(n), parent(parent), d(d), h(h) {}

    Node<N, E>* n;
    Node<N, E>* parent;

    // the cost so far
    C d;

    // the heuristical remaining cost + the cost so far
    C h;

    bool operator<(const RouteNode<N, E, C>& p) const { return h > p.h; }
  };

  template <typename N, typename E, typename C>
  using Settled = std::unordered_map<Node<N, E>*, RouteNode<N, E, C> >;

  template <typename N, typename E, typename C>
  using PQ = std::priority_queue<RouteNode<N, E, C> >;

  template <typename N, typename E, typename C>
  struct CostFunc : public util::graph::CostFunc<N, E, C> {
    virtual ~CostFunc() = default; C operator()(const Edge<N, E>* from, const Node<N, E>* n,
                 const Edge<N, E>* to) const {
      UNUSED(from);
      UNUSED(n);
      UNUSED(to);
      return C();
    };
  };

  template <typename N, typename E, typename C>
  struct HeurFunc : public util::graph::HeurFunc<N, E, C> {
    virtual ~HeurFunc() = default;
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
      const util::graph::CostFunc<N, E, C>& costFunc,
      const util::graph::HeurFunc<N, E, C>&,
      std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Node<N, E>*, NList<N, E>*> resNode);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(const std::set<Node<N, E>*> from,
                            const std::set<Node<N, E>*>& to,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static void relax(RouteNode<N, E, C>& cur, const std::set<Node<N, E>*>& to,
                    const util::graph::CostFunc<N, E, C>& costFunc,
                    const util::graph::HeurFunc<N, E, C>& heurFunc,
                    PQ<N, E, C>& pq);

  template <typename N, typename E, typename C>
  static C relaxFwd(RouteNode<N, E, C>& cur, const std::set<Node<N, E>*>& to,
                    const util::graph::CostFunc<N, E, C>& costFunc,
                    const util::graph::HeurFunc<N, E, C>& heurFunc,
                    PQ<N, E, C>& pq, const Settled<N, E, C>& settledBwd);

  template <typename N, typename E, typename C>
  static C relaxBwd(const std::set<Node<N, E>*>& froms, RouteNode<N, E, C>& cur,
                    const util::graph::CostFunc<N, E, C>& costFunc,
                    const util::graph::HeurFunc<N, E, C>& heurFunc,
                    PQ<N, E, C>& pq, const Settled<N, E, C>& settledFwd);

  template <typename N, typename E, typename C>
  static void buildPath(Node<N, E>* curN, Settled<N, E, C>& settledFwd,
                        Settled<N, E, C>& settledBwd, NList<N, E>* resNodes,
                        EList<N, E>* resEdges);

  static size_t ITERS;
};

#include "util/graph/BiDijkstra.tpp"
}  // namespace graph
}  // namespace util

#endif  // UTIL_GRAPH_BIDIJKSTRA_H_
