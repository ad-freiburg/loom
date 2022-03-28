// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_EDIJKSTRA_H_
#define UTIL_GRAPH_EDIJKSTRA_H_

#include <limits>
#include <list>
#include <set>
#include <unordered_map>
#include "util/PriorityQueue.h"
#include "util/graph/Edge.h"
#include "util/graph/Graph.h"
#include "util/graph/Node.h"
#include "util/graph/ShortestPath.h"
#include "util/graph/radix_heap.h"
#include "util/graph/robin/robin_map.h"

namespace util {
namespace graph {

using util::graph::Edge;
using util::graph::Graph;
using util::graph::Node;

// edge-based dijkstra - settles edges instead of nodes
class EDijkstra : public ShortestPath<EDijkstra> {
 public:
  template <typename N, typename E, typename C>
  struct RouteEdge {
    RouteEdge() : e(0), parent(0), d(), n(0) {}
    RouteEdge(Edge<N, E>* e) : e(e), parent(0), d(), n(0) {}
    RouteEdge(Edge<N, E>* e, Edge<N, E>* parent, Node<N, E>* n, C d)
        : e(e), parent(parent), d(d), n(n) {}

    Edge<N, E>* e;
    Edge<N, E>* parent;

    // the cost so far
    C d;

    Node<N, E>* n;
  };

  template <typename N, typename E, typename C>
  struct RouteEdgeInit {
    RouteEdgeInit() : e(0), parent(0), d(), dwi(), n(0) {}
    RouteEdgeInit(Edge<N, E>* e) : e(e), parent(0), d(), dwi(), n(0) {}
    RouteEdgeInit(Edge<N, E>* e, Edge<N, E>* parent, Node<N, E>* n, C d)
        : e(e), parent(parent), d(d), dwi(), n(n) {}
    RouteEdgeInit(Edge<N, E>* e, Edge<N, E>* parent, Node<N, E>* n, C d, C dwi)
        : e(e), parent(parent), d(d), dwi(dwi), n(n) {}

    Edge<N, E>* e;
    Edge<N, E>* parent;

    // the cost so far
    C d;

    // the cost without the initial costs
    C dwi;

    Node<N, E>* n;
  };

  template <typename N, typename E, typename C>
  struct RouteEdgeInitNoRes {
    RouteEdgeInitNoRes() : e(0), parent(0), d(), dwi() {}
    RouteEdgeInitNoRes(Edge<N, E>* e) : e(e), parent(0), d(), dwi() {}
    RouteEdgeInitNoRes(Edge<N, E>* e, Edge<N, E>* parent, C d)
        : e(e), parent(parent), d(d), dwi() {}
    RouteEdgeInitNoRes(Edge<N, E>* e, Edge<N, E>* parent, C d, C dwi)
        : e(e), parent(parent), d(d), dwi(dwi) {}

    Edge<N, E>* e;
    Edge<N, E>* parent;

    // the cost so far
    C d;

    // the cost without the initial costs
    C dwi;
  };

  template <typename N, typename E, typename C>
  struct CostFunc : public util::graph::CostFunc<N, E, C> {
    C operator()(const Node<N, E>* from, const Edge<N, E>* e,
                 const Node<N, E>* to) const {
      UNUSED(from);
      UNUSED(e);
      UNUSED(to);
      return C();
    };
  };

  template <typename N, typename E, typename C>
  struct HeurFunc : public util::graph::HeurFunc<N, E, C> {
    C operator()(const Node<N, E>* from,
                 const std::set<Node<N, E>*>& to) const {
      UNUSED(from);
      UNUSED(to);
      return C();
    };
  };

  template <typename N, typename E, typename C>
  using Settled = tsl::robin_map<Edge<N, E>*, RouteEdge<N, E, C>>;

  template <typename N, typename E, typename C>
  using PQ = radix_heap::pair_radix_heap<C, RouteEdge<N, E, C>>;

  template <typename N, typename E, typename C>
  using SettledInit = tsl::robin_map<Edge<N, E>*, RouteEdgeInit<N, E, C>>;

  template <typename N, typename E, typename C>
  using SettledInitNoRes =
      tsl::robin_map<Edge<N, E>*, RouteEdgeInitNoRes<N, E, C>>;

  template <typename N, typename E, typename C>
  using PQInit = radix_heap::pair_radix_heap<C, RouteEdgeInit<N, E, C>>;

  template <typename N, typename E, typename C>
  using PQInitNoRes =
      radix_heap::pair_radix_heap<C, RouteEdgeInitNoRes<N, E, C>>;

  template <typename N, typename E, typename C>
  static C shortestPathImpl(const std::set<Edge<N, E>*>& from,
                            const std::set<Edge<N, E>*>& to,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(Edge<N, E>* from, const std::set<Node<N, E>*>& to,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(const std::set<Edge<N, E>*>& from,
                            const std::set<Node<N, E>*>& to,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static C shortestPathImpl(const std::set<Node<N, E>*>& from,
                            const std::set<Node<N, E>*>& to,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc,
                            EList<N, E>* resEdges, NList<N, E>* resNodes);

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathImpl(
      const std::set<Edge<N, E>*>& from,
      const util::graph::CostFunc<N, E, C>& costFunc, bool rev);

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathImpl(
      Edge<N, E>* from, const std::set<Edge<N, E>*>& to,
      const util::graph::CostFunc<N, E, C>& costFunc,
      const util::graph::HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes);

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathImpl(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts,
      const util::graph::CostFunc<N, E, C>& costFunc,
      const util::graph::HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes);

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathImpl(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts, C stall,
      const util::graph::CostFunc<N, E, C>& costFunc,
      const util::graph::HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes);

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, std::pair<Edge<N, E>*, C>>
  shortestPathImpl(const std::set<Edge<N, E>*>& from,
                   const std::set<Edge<N, E>*>& to,
                   const std::unordered_map<Edge<N, E>*, C>& initCosts, C stall,
                   const util::graph::CostFunc<N, E, C>& costFunc,
                   const util::graph::HeurFunc<N, E, C>& heurFunc);

  template <typename N, typename E, typename C>
  static void buildPath(Edge<N, E>* curE, const Settled<N, E, C>& settled,
                        NList<N, E>* resNodes, EList<N, E>* resEdges);

  template <typename N, typename E, typename C>
  static void buildPathInit(Edge<N, E>* curE,
                            const SettledInit<N, E, C>& settled,
                            NList<N, E>* resNodes, EList<N, E>* resEdges);

  template <typename N, typename E, typename C>
  static inline void relax(RouteEdge<N, E, C>& cur,
                           const std::set<Edge<N, E>*>& to,
                           const util::graph::CostFunc<N, E, C>& costFunc,
                           const util::graph::HeurFunc<N, E, C>& heurFunc,
                           PQ<N, E, C>& pq);

  template <typename N, typename E, typename C>
  static inline void relaxInit(RouteEdgeInit<N, E, C>& cur,
                               const std::set<Edge<N, E>*>& to, C stall,
                               const util::graph::CostFunc<N, E, C>& costFunc,
                               const util::graph::HeurFunc<N, E, C>& heurFunc,
                               PQInit<N, E, C>& pq);

  template <typename N, typename E, typename C>
  static inline void relaxInitNoResEdgs(
      RouteEdgeInitNoRes<N, E, C>& cur, const std::set<Edge<N, E>*>& to,
      C stall, const util::graph::CostFunc<N, E, C>& costFunc,
      const util::graph::HeurFunc<N, E, C>& heurFunc, PQInitNoRes<N, E, C>& pq);

  template <typename N, typename E, typename C>
  static void relaxInv(RouteEdge<N, E, C>& cur,
                       const util::graph::CostFunc<N, E, C>& costFunc,
                       PQ<N, E, C>& pq);

  static size_t ITERS;
};

#include "util/graph/EDijkstra.tpp"
}  // namespace graph
}  // namespace util

#endif  // UTIL_GRAPH_DIJKSTRA_H_
