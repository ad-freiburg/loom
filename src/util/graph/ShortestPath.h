// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GRAPH_SHORTESTPATH_H_
#define UTIL_GRAPH_SHORTESTPATH_H_

#include <exception>
#include <iostream>
#include <limits>
#include <list>
#include <queue>
#include <set>
#include <unordered_map>
#include "util/graph/Edge.h"
#include "util/graph/Graph.h"
#include "util/graph/Node.h"

namespace util {
namespace graph {

using util::graph::Edge;
using util::graph::Graph;
using util::graph::Node;

template <typename N, typename E>
using EList = std::vector<Edge<N, E>*>;

template <typename N, typename E>
using NList = std::vector<Node<N, E>*>;

template <typename N, typename E, typename C>
struct CostFunc {
  virtual C operator()(const Node<N, E>* from, const Edge<N, E>* e,
                       const Node<N, E>* to) const = 0;
  virtual C operator()(const Edge<N, E>* from, const Node<N, E>* n,
                       const Edge<N, E>* to) const = 0;
  virtual C inf() const = 0;
};

template <typename N, typename E, typename C>
struct HeurFunc {
  virtual C operator()(const Node<N, E>* a,
                       const std::set<Node<N, E>*>& b) const = 0;
  virtual C operator()(const Edge<N, E>* a,
                       const std::set<Edge<N, E>*>& b) const = 0;
};

template <typename N, typename E, typename C>
struct ZeroHeurFunc : public HeurFunc<N, E, C> {
  C operator()(const Node<N, E>* a, const std::set<Node<N, E>*>& b) const {
    UNUSED(a);
    UNUSED(b);
    return C();
  }
  C operator()(const Edge<N, E>* a, const std::set<Edge<N, E>*>& b) const {
    UNUSED(a);
    UNUSED(b);
    return C();
  }
};

// shortest path base class
template <class D>
class ShortestPath {
 public:
  template <typename N, typename E>
  using EList = std::vector<Edge<N, E>*>;

  template <typename N, typename E>
  using NList = std::vector<Node<N, E>*>;

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        EList<N, E>* resEdges, NList<N, E>* resNodes) {
    return D::shortestPathImpl(from, to, costFunc, heurFunc, resEdges,
                               resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(const std::set<Node<N, E>*>& from,
                        const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc) {
    EList<N, E>* el = 0;
    NList<N, E>* nl = 0;
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(),
                               el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(const std::set<Node<N, E>*> from,
                        const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        EList<N, E>* resEdges, NList<N, E>* resNodes) {
    return D::shortestPathImpl(from, to, costFunc, heurFunc, resEdges,
                               resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(const std::set<Node<N, E>*> from,
                        const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        EList<N, E>* resEdges, NList<N, E>* resNodes) {
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(),
                               resEdges, resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        EList<N, E>* resEdges, NList<N, E>* resNodes) {
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(),
                               resEdges, resNodes);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Node<N, E>*, C> shortestPath(
      Node<N, E>* from, const std::set<Node<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Node<N, E>*, NList<N, E>*> resNodes) {
    return D::shortestPathImpl(from, to, costFunc, heurFunc, resEdges,
                               resNodes);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Node<N, E>*, C> shortestPath(
      Node<N, E>* from, const std::set<Node<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc,
      std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Node<N, E>*, NList<N, E>*> resNodes) {
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(),
                               resEdges, resNodes);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Node<N, E>*, C> shortestPath(
      Node<N, E>* from, const std::set<Node<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges) {
    std::unordered_map<Node<N, E>*, NList<N, E>*> dummyRet;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, resEdges,
                               dummyRet);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Node<N, E>*, C> shortestPath(
      Node<N, E>* from, const std::set<Node<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc,
      std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges) {
    std::unordered_map<Node<N, E>*, NList<N, E>*> dummyRet;
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(),
                               resEdges, dummyRet);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Node<N, E>*, C> shortestPath(
      Node<N, E>* from, const std::set<Node<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc,
      std::unordered_map<Node<N, E>*, NList<N, E>*> resNodes) {
    std::unordered_map<Node<N, E>*, EList<N, E>*> dummyRet;
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(),
                               dummyRet, resNodes);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Node<N, E>*, C> shortestPath(
      Node<N, E>* from, const std::set<Node<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Node<N, E>*, NList<N, E>*> resNodes) {
    std::unordered_map<Node<N, E>*, EList<N, E>*> dummyRet;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, dummyRet,
                               resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        EList<N, E>* resEdges, NList<N, E>* resNodes) {
    if (to->getInDeg() == 0) return costFunc.inf();
    std::set<Node<N, E>*> tos;
    tos.insert(to);
    return shortestPath(from, tos, costFunc, resEdges, resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, const std::set<Node<N, E>*>& tos,
                        const CostFunc<N, E, C>& costFunc) {
    EList<N, E>* el = 0;
    NList<N, E>* nl = 0;
    return shortestPath(from, tos, costFunc, el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc) {
    if (to->getInDeg() == 0) return costFunc.inf();
    std::set<Node<N, E>*> tos;
    tos.insert(to);
    EList<N, E>* el = 0;
    NList<N, E>* nl = 0;
    return shortestPath(from, tos, costFunc, el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        NList<N, E>* resNodes) {
    if (to->getInDeg() == 0) return costFunc.inf();
    std::set<Node<N, E>*> tos;
    tos.insert(to);
    EList<N, E>* el = 0;
    return shortestPath(from, tos, costFunc, el, resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        EList<N, E>* resEdges) {
    if (to->getInDeg() == 0) return costFunc.inf();
    std::set<Node<N, E>*> tos;
    tos.insert(to);
    NList<N, E>* nl = 0;
    return shortestPath(from, tos, costFunc, resEdges, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        NList<N, E>* resNodes) {
    EList<N, E>* el = 0;
    return shortestPath(from, to, costFunc, el, resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        EList<N, E>* resEdges) {
    NList<N, E>* nl = 0;
    return shortestPath(from, to, costFunc, resEdges, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        EList<N, E>* resEdges, NList<N, E>* resNodes) {
    if (to->getInDeg() == 0) return costFunc.inf();
    std::set<Node<N, E>*> tos;
    tos.insert(to);
    return D::shortestPathImpl(from, tos, costFunc, heurFunc, resEdges,
                               resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        NList<N, E>* resNodes) {
    if (to->getInDeg() == 0) return costFunc.inf();
    std::set<Node<N, E>*> tos;
    tos.insert(to);
    EList<N, E>* el = 0;
    return D::shortestPathImpl(from, tos, costFunc, heurFunc, el, resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        EList<N, E>* resEdges) {
    if (to->getInDeg() == 0) return costFunc.inf();
    std::set<Node<N, E>*> tos;
    tos.insert(to);
    NList<N, E>* nl = 0;
    return D::shortestPathImpl(from, tos, costFunc, heurFunc, resEdges, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        NList<N, E>* resNodes) {
    EList<N, E>* el = 0;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, el, resNodes);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        EList<N, E>* resEdges) {
    NList<N, E>* nl = 0;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, resEdges, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from, Edge<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc,
                        EList<N, E>* resEdges) {
    NList<N, E> dummyRet;
    std::set<Edge<N, E>*> froms{from};
    std::set<Edge<N, E>*> tos{to};
    return D::shortestPathImpl(froms, tos, costFunc, heurFunc, resEdges,
                               &dummyRet);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      Edge<N, E>* from, const std::set<Edge<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes) {
    return D::shortestPathImpl(from, to, costFunc, heurFunc, resEdges,
                               resNodes);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      Edge<N, E>* from, const std::set<Edge<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes) {
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>()
, resEdges,
                               resNodes);
  }



  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      Edge<N, E>* from, const std::set<Edge<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges) {
    std::unordered_map<Edge<N, E>*, NList<N, E>*> dummyRet;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, resEdges,
                               dummyRet);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      Edge<N, E>* from, const std::set<Edge<N, E>*>& to,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc) {
    std::unordered_map<Edge<N, E>*, NList<N, E>*> dummyRet;
    std::unordered_map<Edge<N, E>*, EList<N, E>*> dummyRetE;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, dummyRetE,
                               dummyRet);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from,
                        const std::set<Edge<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc) {
    NList<N, E>* nl = 0;
    EList<N, E>* el = 0;
    std::set<Edge<N, E>*> fromS;
    fromS.insert(from);
    return D::shortestPathImpl(fromS, to, costFunc, ZeroHeurFunc<N, E, C>(), el,
                               nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(const std::set<Edge<N, E>*>& from,
                        const std::set<Edge<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc) {
    NList<N, E>* nl = 0;
    EList<N, E>* el = 0;
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(), el,
                               nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(const std::set<Edge<N, E>*>& from,
                        const std::set<Edge<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc) {
    NList<N, E>* nl = 0;
    EList<N, E>* el = 0;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(const std::set<Edge<N, E>*>& from,
                        const std::set<Edge<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc, EList<N, E>* el) {
    NList<N, E>* nl = 0;
    return D::shortestPathImpl(from, to, costFunc, heurFunc, el, nl);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes) {
    return D::shortestPathImpl(from, to, initCosts, costFunc.inf(), costFunc,
                               heurFunc, resEdges, resNodes);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges) {
    std::unordered_map<Edge<N, E>*, NList<N, E>*> dummyRet;
    return D::shortestPathImpl(from, to, initCosts, costFunc.inf(), costFunc,
                               heurFunc, resEdges, dummyRet);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc) {
    std::unordered_map<Edge<N, E>*, NList<N, E>*> dummyRet;
    std::unordered_map<Edge<N, E>*, EList<N, E>*> dummyRetE;
    return D::shortestPathImpl(from, to, initCosts, costFunc.inf(), costFunc,
                               heurFunc, dummyRetE, dummyRet);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts, C stall,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
      std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes) {
    return D::shortestPathImpl(from, to, initCosts, stall, costFunc, heurFunc,
                               resEdges, resNodes);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts, C stall,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc,
      std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges) {
    std::unordered_map<Edge<N, E>*, NList<N, E>*> dummyRet;
    return D::shortestPathImpl(from, to, initCosts, stall, costFunc, heurFunc,
                               resEdges, dummyRet);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, std::pair<Edge<N, E>*, C>> shortestPath(
      const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
      const std::unordered_map<Edge<N, E>*, C>& initCosts, C stall,
      const CostFunc<N, E, C>& costFunc, const HeurFunc<N, E, C>& heurFunc) {
    return D::shortestPathImpl(from, to, initCosts, stall, costFunc, heurFunc);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from, Edge<N, E>* to,
                        const CostFunc<N, E, C>& costFunc) {
    NList<N, E>* nl = 0;
    EList<N, E>* el = 0;
    std::set<Edge<N, E>*> tos{to};
    std::set<Edge<N, E>*> froms{from};
    return D::shortestPathImpl(froms, tos, costFunc, ZeroHeurFunc<N, E, C>(),
                               el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from, Edge<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc) {
    NList<N, E>* nl = 0;
    EList<N, E>* el = 0;
    std::set<Edge<N, E>*> tos{to};
    std::set<Edge<N, E>*> froms{from};
    return D::shortestPathImpl(froms, tos, costFunc, heurFunc, el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc, EList<N, E>* el,
                        NList<N, E>* nl) {
    return D::shortestPathImpl(from, to, costFunc, heurFunc, el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from, const std::set<Node<N, E>*>& to,
                        const CostFunc<N, E, C>& costFunc, EList<N, E>* el,
                        NList<N, E>* nl) {
    return D::shortestPathImpl(from, to, costFunc, ZeroHeurFunc<N, E, C>(), el,
                               nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc,
                        const HeurFunc<N, E, C>& heurFunc, EList<N, E>* el,
                        NList<N, E>* nl) {
    std::set<Node<N, E>*> tos{to};
    return D::shortestPathImpl(from, tos, costFunc, heurFunc, el, nl);
  }

  template <typename N, typename E, typename C>
  static C shortestPath(Edge<N, E>* from, Node<N, E>* to,
                        const CostFunc<N, E, C>& costFunc, EList<N, E>* el,
                        NList<N, E>* nl) {
    std::set<Node<N, E>*> tos{to};
    return D::shortestPathImpl(from, tos, costFunc, ZeroHeurFunc<N, E, C>(), el,
                               nl);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      Edge<N, E>* from, const CostFunc<N, E, C>& costFunc) {
    std::set<Edge<N, E>*> froms{from};
    return D::shortestPathImpl(froms, costFunc, false);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathRev(
      Edge<N, E>* from, const CostFunc<N, E, C>& costFunc) {
    std::set<Edge<N, E>*> froms{from};
    return D::shortestPathImpl(froms, costFunc, true);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPath(
      Node<N, E>* from, const CostFunc<N, E, C>& costFunc) {
    std::set<Edge<N, E>*> froms;
    froms.insert(from->getAdjListOut().begin(), from->getAdjListOut().end());
    return D::shortestPathImpl(froms, costFunc, false);
  }

  template <typename N, typename E, typename C>
  static std::unordered_map<Edge<N, E>*, C> shortestPathRev(
      Node<N, E>* from, const CostFunc<N, E, C>& costFunc) {
    std::set<Edge<N, E>*> froms;
    froms.insert(from->getAdjListOut().begin(), from->getAdjListOut().end());
    return D::shortestPathImpl(froms, costFunc, true);
  }
};
}  // namespace graph
}  // namespace util

#endif  // UTIL_GRAPH_SHORTESTPATH_H_
