// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C EDijkstra::shortestPathImpl(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                              const util::graph::CostFunc<N, E, C>& costFunc,
                              const util::graph::HeurFunc<N, E, C>& heurFunc,
                              EList<N, E>* resEdges, NList<N, E>* resNodes) {
  std::set<Edge<N, E>*> frEs;
  std::set<Edge<N, E>*> toEs;

  frEs.insert(from->getAdjListOut().begin(), from->getAdjListOut().end());

  for (auto n : to) {
    toEs.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
  }

  C cost = shortestPathImpl(frEs, toEs, costFunc, heurFunc, resEdges, resNodes);

  // the beginning node is not included in our edge based dijkstra
  if (resNodes) resNodes->push_back(from);

  return cost;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C EDijkstra::shortestPathImpl(Edge<N, E>* from, const std::set<Node<N, E>*>& to,
                              const util::graph::CostFunc<N, E, C>& costFunc,
                              const util::graph::HeurFunc<N, E, C>& heurFunc,
                              EList<N, E>* resEdges, NList<N, E>* resNodes) {
  std::set<Edge<N, E>*> frEs;
  std::set<Edge<N, E>*> toEs;

  frEs.insert(from);

  for (auto n : to) {
    toEs.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
  }

  C cost = shortestPathImpl(frEs, toEs, costFunc, heurFunc, resEdges, resNodes);

  return cost;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C EDijkstra::shortestPathImpl(const std::set<Node<N, E>*>& from,
                              const std::set<Node<N, E>*>& to,
                              const util::graph::CostFunc<N, E, C>& costFunc,
                              const util::graph::HeurFunc<N, E, C>& heurFunc,
                              EList<N, E>* resEdges, NList<N, E>* resNodes) {
  std::set<Edge<N, E>*> frEs;
  std::set<Edge<N, E>*> toEs;

  for (auto n : from) {
    frEs.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
  }

  for (auto n : to) {
    toEs.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
  }

  C cost = shortestPathImpl(frEs, toEs, costFunc, heurFunc, resEdges, resNodes);

  return cost;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C EDijkstra::shortestPathImpl(const std::set<Edge<N, E>*>& from,
                              const std::set<Edge<N, E>*>& to,
                              const util::graph::CostFunc<N, E, C>& costFunc,
                              const util::graph::HeurFunc<N, E, C>& heurFunc,
                              EList<N, E>* resEdges, NList<N, E>* resNodes) {
  if (from.size() == 0 || to.size() == 0) return costFunc.inf();

  Settled<N, E, C> settled;
  PQ<N, E, C> pq;
  bool found = false;

  // at the beginning, put all edges on the priority queue,
  // init them with their own cost
  for (auto e : from) {
    C c = costFunc(0, 0, e);
    C h = heurFunc(e, to);
    pq.push(c + h, {e, (Edge<N, E>*)0, (Node<N, E>*)0, c});
  }

  RouteEdge<N, E, C> cur;

  while (!pq.empty()) {
    if (costFunc.inf() <= pq.topKey()) return costFunc.inf();
    auto se = settled.find(pq.topVal().e);
    if (se != settled.end()) {
      // to allow non-consistent heuristics
      if (se->second.d <= pq.topVal().d) {
        pq.pop();
        continue;
      }
    }
    EDijkstra::ITERS++;

    cur = pq.topVal();
    pq.pop();

    settled[cur.e] = cur;

    if (to.find(cur.e) != to.end()) {
      found = true;
      break;
    }

    relax(cur, to, costFunc, heurFunc, pq);
  }

  if (!found) return costFunc.inf();

  buildPath(cur.e, settled, resNodes, resEdges);

  return cur.d;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
std::unordered_map<Edge<N, E>*, C> EDijkstra::shortestPathImpl(
    const std::set<Edge<N, E>*>& from,
    const util::graph::CostFunc<N, E, C>& costFunc, bool rev) {
  std::unordered_map<Edge<N, E>*, C> costs;

  Settled<N, E, C> settled;
  PQ<N, E, C> pq;

  std::set<Edge<N, E>*> to;

  for (auto e : from) {
    pq.push(C(), {e, (Edge<N, E>*)0, (Node<N, E>*)0, costFunc(0, 0, e)});
  }

  RouteEdge<N, E, C> cur;

  while (!pq.empty()) {
    auto se = settled.find(pq.topVal().e);
    if (se != settled.end()) {
      // to allow non-consistent heuristics
      if (se->second.d <= pq.topVal().d) {
        pq.pop();
        continue;
      }
    }
    EDijkstra::ITERS++;

    cur = pq.topVal();
    pq.pop();

    settled[cur.e] = cur;

    costs[cur.e] = cur.d;
    buildPath(cur.e, settled, (NList<N, E>*)0, (EList<N, E>*)0);

    if (rev)
      relaxInv(cur, costFunc, pq);
    else
      relax(cur, to, costFunc, ZeroHeurFunc<N, E, C>(), pq);
  }

  return costs;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
std::unordered_map<Edge<N, E>*, C> EDijkstra::shortestPathImpl(
    Edge<N, E>* from, const std::set<Edge<N, E>*>& to,
    const util::graph::CostFunc<N, E, C>& costFunc,
    const util::graph::HeurFunc<N, E, C>& heurFunc,
    std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
    std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes) {
  std::unordered_map<Edge<N, E>*, C> costs;
  if (to.size() == 0) return costs;

  // init costs with inf
  for (auto e : to) costs[e] = costFunc.inf();

  Settled<N, E, C> settled;
  PQ<N, E, C> pq;

  size_t found = 0;

  C c = costFunc(0, 0, from);
  C h = heurFunc(from, to);
  pq.push(c + h, {from, (Edge<N, E>*)0, (Node<N, E>*)0, c});

  RouteEdge<N, E, C> cur;

  while (!pq.empty()) {
    if (costFunc.inf() <= pq.topKey()) return costs;
    auto se = settled.find(pq.topVal().e);
    if (se != settled.end()) {
      // to allow non-consistent heuristics
      if (se->second.d <= pq.topVal().d) {
        pq.pop();
        continue;
      }
    }
    EDijkstra::ITERS++;

    cur = pq.topVal();
    pq.pop();

    settled[cur.e] = cur;

    if (to.find(cur.e) != to.end()) {
      found++;
      costs[cur.e] = cur.d;
      buildPath(cur.e, settled, resNodes[cur.e], resEdges[cur.e]);
    }

    if (found == to.size()) return costs;

    relax(cur, to, costFunc, heurFunc, pq);
  }

  return costs;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
std::unordered_map<Edge<N, E>*, C> EDijkstra::shortestPathImpl(
    const std::set<Edge<N, E>*>& from, const std::set<Edge<N, E>*>& to,
    const std::unordered_map<Edge<N, E>*, C>& initCosts, C stall,
    const util::graph::CostFunc<N, E, C>& costFunc,
    const util::graph::HeurFunc<N, E, C>& heurFunc,
    std::unordered_map<Edge<N, E>*, EList<N, E>*> resEdges,
    std::unordered_map<Edge<N, E>*, NList<N, E>*> resNodes) {
  /**
   * Shortest paths from the set <from> to ALL nodes in <TO>, but
   * init <from> nodes with costs (this is equivalent to adding an auxiliary
   * node S, connecting it with directed edges to all <from>, setting the
   * costs of these edges to the initial costs and run a 1->N Dijkstra from S
   **/

  std::unordered_map<Edge<N, E>*, C> costs;
  if (to.size() == 0) return costs;

  // init costs with inf
  for (auto e : to) costs[e] = costFunc.inf();

  SettledInit<N, E, C> settled;
  PQInit<N, E, C> pq;

  size_t found = 0;

  // put all nodes in from onto the PQ with their initial costs, also set
  // the initial cost as a heuristic starting point!
  for (auto e : from) {
    C iCost = initCosts.find(e)->second;
    assert(iCost + heurFunc(e, to) >= iCost);
    pq.push(iCost + heurFunc(e, to),
            {e, (Edge<N, E>*)0, (Node<N, E>*)0, iCost});
  }

  RouteEdgeInit<N, E, C> cur;

  while (!pq.empty()) {
    if (costFunc.inf() <= pq.topKey()) return costs;
    if (stall <= pq.topVal().dwi) return costs;
    auto se = settled.find(pq.topVal().e);
    if (se != settled.end()) {
      // to allow non-consistent heuristics
      if (se->second.d <= pq.topVal().d) {
        pq.pop();
        continue;
      }
    }
    EDijkstra::ITERS++;

    cur = pq.topVal();
    pq.pop();

    settled[cur.e] = cur;

    if (to.find(cur.e) != to.end()) {
      found++;
      costs[cur.e] = cur.d;
      buildPathInit(cur.e, settled, resNodes[cur.e], resEdges[cur.e]);

      if (found == to.size()) return costs;
    }

    relaxInit(cur, to, stall, costFunc, heurFunc, pq);
  }

  return costs;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
std::unordered_map<Edge<N, E>*, std::pair<Edge<N, E>*, C>>
EDijkstra::shortestPathImpl(const std::set<Edge<N, E>*>& from,
                            const std::set<Edge<N, E>*>& to,
                            const std::unordered_map<Edge<N, E>*, C>& initCosts,
                            C stall,
                            const util::graph::CostFunc<N, E, C>& costFunc,
                            const util::graph::HeurFunc<N, E, C>& heurFunc) {
  /**
   * Shortest paths from the set <from> to ALL nodes in <TO>, but
   * init <from> nodes with costs (this is equivalent to adding an auxiliary
   * node S, connecting it with directed edges to all <from>, setting the
   * costs of these edges to the initial costs and run a 1->N Dijkstra from S
   **/

  std::unordered_map<Edge<N, E>*, std::pair<Edge<N, E>*, C>> costs;
  if (to.size() == 0) return costs;

  // init costs with inf
  for (auto e : to) costs[e] = {0, costFunc.inf()};

  SettledInitNoRes<N, E, C> settled;
  PQInitNoRes<N, E, C> pq;

  size_t found = 0;

  // put all nodes in from onto the PQ with their initial costs, also set
  // the initial cost as a heuristic starting point!
  // set the parent to the edge itself - in this version, the parent is ALWAYS
  // the start edge, as we don't need the exact paths later on
  for (auto e : from) {
    C iCost = initCosts.find(e)->second;
    assert(iCost + heurFunc(e, to) >= iCost);
    pq.push(iCost + heurFunc(e, to), {e, e, iCost});
  }

  RouteEdgeInitNoRes<N, E, C> cur;

  while (!pq.empty()) {
    if (costFunc.inf() <= pq.topKey()) return costs;
    if (stall <= pq.topVal().dwi) return costs;
    auto se = settled.find(pq.topVal().e);
    if (se != settled.end()) {
      // to allow non-consistent heuristics
      if (se->second.d <= pq.topVal().d) {
        pq.pop();
        continue;
      }
    }
    EDijkstra::ITERS++;

    cur = pq.topVal();
    pq.pop();

    settled[cur.e] = cur;

    if (to.find(cur.e) != to.end()) {
      found++;
      costs[cur.e] = {cur.parent, cur.d};
      if (found == to.size()) return costs;
    }

    relaxInitNoResEdgs(cur, to, stall, costFunc, heurFunc, pq);
  }

  return costs;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::relaxInv(RouteEdge<N, E, C>& cur,
                         const util::graph::CostFunc<N, E, C>& costFunc,
                         PQ<N, E, C>& pq) {
  // handling undirected graph makes no sense here

  for (const auto edge : cur.e->getFrom()->getAdjListIn()) {
    if (edge == cur.e) continue;
    C newC = costFunc(edge, cur.e->getFrom(), cur.e);
    newC = cur.d + newC;
    if (costFunc.inf() <= newC) continue;
    if (newC < cur.d) continue;  // cost overflow!

    pq.push(C(), {edge, cur.e, cur.e->getFrom(), newC});
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::relaxInitNoResEdgs(
    RouteEdgeInitNoRes<N, E, C>& cur, const std::set<Edge<N, E>*>& to, C stall,
    const util::graph::CostFunc<N, E, C>& costFunc,
    const util::graph::HeurFunc<N, E, C>& heurFunc, PQInitNoRes<N, E, C>& pq) {
  if (cur.e->getFrom()->hasEdgeIn(cur.e) &&
      cur.e->getFrom() != cur.e->getTo()) {
    // for undirected graphs
    for (const auto edge : cur.e->getFrom()->getAdjListOut()) {
      if (edge == cur.e) continue;
      C newC = costFunc(cur.e, cur.e->getFrom(), edge);
      C newDwi = cur.dwi + newC;

      if (stall <= newDwi) continue;

      if (costFunc.inf() <= newC) continue;
      newC = cur.d + newC;
      if (newC < cur.d) continue;  // cost overflow!

      if (costFunc.inf() <= newC) continue;

      const C& h = heurFunc(edge, to);
      if (costFunc.inf() <= h) continue;
      newC = cur.d + newC;
      const C& newH = newC + h;
      if (newH < newC) continue;  // cost overflow!

      pq.push(newH, {edge, cur.parent, newC, newDwi});
    }
  }

  for (const auto edge : cur.e->getTo()->getAdjListOut()) {
    if (edge == cur.e) continue;
    C newC = costFunc(cur.e, cur.e->getTo(), edge);
    C newDwi = cur.dwi + newC;

    if (stall <= newDwi) continue;

    if (costFunc.inf() <= newC) continue;
    newC = cur.d + newC;
    if (newC < cur.d) continue;  // cost overflow!

    if (costFunc.inf() <= newC) continue;

    const C& h = heurFunc(edge, to);
    if (costFunc.inf() <= h) continue;
    const C& newH = newC + h;
    if (newH < newC) continue;  // cost overflow!

    pq.push(newH, {edge, cur.parent, newC, newDwi});
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::relaxInit(RouteEdgeInit<N, E, C>& cur,
                          const std::set<Edge<N, E>*>& to, C stall,
                          const util::graph::CostFunc<N, E, C>& costFunc,
                          const util::graph::HeurFunc<N, E, C>& heurFunc,
                          PQInit<N, E, C>& pq) {
  if (cur.e->getFrom()->hasEdgeIn(cur.e) &&
      cur.e->getFrom() != cur.e->getTo()) {
    // for undirected graphs
    for (const auto edge : cur.e->getFrom()->getAdjListOut()) {
      if (edge == cur.e) continue;
      C newC = costFunc(cur.e, cur.e->getFrom(), edge);
      C newDwi = cur.dwi + newC;

      if (stall <= newDwi) continue;

      if (costFunc.inf() <= newC) continue;
      newC = cur.d + newC;
      if (newC < cur.d) continue;  // cost overflow!

      if (costFunc.inf() <= newC) continue;

      const C& h = heurFunc(edge, to);
      if (costFunc.inf() <= h) continue;
      newC = cur.d + newC;
      const C& newH = newC + h;
      if (newH < newC) continue;  // cost overflow!

      pq.push(newH, {edge, cur.e, cur.e->getFrom(), newC, newDwi});
    }
  }

  for (const auto edge : cur.e->getTo()->getAdjListOut()) {
    if (edge == cur.e) continue;
    C newC = costFunc(cur.e, cur.e->getTo(), edge);
    C newDwi = cur.dwi + newC;

    if (stall <= newDwi) continue;

    if (costFunc.inf() <= newC) continue;
    newC = cur.d + newC;
    if (newC < cur.d) continue;  // cost overflow!

    if (costFunc.inf() <= newC) continue;

    const C& h = heurFunc(edge, to);
    if (costFunc.inf() <= h) continue;
    const C& newH = newC + h;
    if (newH < newC) continue;  // cost overflow!

    pq.push(newH, {edge, cur.e, cur.e->getTo(), newC, newDwi});
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::relax(RouteEdge<N, E, C>& cur, const std::set<Edge<N, E>*>& to,
                      const util::graph::CostFunc<N, E, C>& costFunc,
                      const util::graph::HeurFunc<N, E, C>& heurFunc,
                      PQ<N, E, C>& pq) {
  if (cur.e->getFrom()->hasEdgeIn(cur.e) &&
      cur.e->getFrom() != cur.e->getTo()) {
    // for undirected graphs
    for (const auto edge : cur.e->getFrom()->getAdjListOut()) {
      if (edge == cur.e) continue;
      C newC = costFunc(cur.e, cur.e->getFrom(), edge);
      if (costFunc.inf() <= newC) continue;
      newC = cur.d + newC;
      if (newC < cur.d) continue;  // cost overflow!
      if (costFunc.inf() <= newC) continue;

      const C& h = heurFunc(edge, to);
      if (costFunc.inf() <= h) continue;
      const C& newH = newC + h;
      if (newH < newC) continue;  // cost overflow!

      pq.push(newH, {edge, cur.e, cur.e->getFrom(), newC});
    }
  }

  for (const auto edge : cur.e->getTo()->getAdjListOut()) {
    if (edge == cur.e) continue;
    C newC = costFunc(cur.e, cur.e->getTo(), edge);
    if (costFunc.inf() <= newC) continue;
    newC = cur.d + newC;
    if (newC < cur.d) continue;  // cost overflow!
    if (costFunc.inf() <= newC) continue;

    const C& h = heurFunc(edge, to);
    if (costFunc.inf() <= h) continue;
    const C& newH = newC + h;
    if (newH < newC) continue;  // cost overflow!

    pq.push(newH, {edge, cur.e, cur.e->getTo(), newC});
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::buildPath(Edge<N, E>* curE, const Settled<N, E, C>& settled,
                          NList<N, E>* resNodes, EList<N, E>* resEdges) {
  const RouteEdge<N, E, C>* curEdge = &settled.find(curE)->second;
  if (resNodes) resNodes->push_back(curEdge->e->getOtherNd(curEdge->n));
  while (resNodes || resEdges) {
    if (resNodes && curEdge->n) resNodes->push_back(curEdge->n);
    if (resEdges) resEdges->push_back(curEdge->e);
    if (!curEdge->parent) break;
    curEdge = &settled.find(curEdge->parent)->second;
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::buildPathInit(Edge<N, E>* curE,
                              const SettledInit<N, E, C>& settled,
                              NList<N, E>* resNodes, EList<N, E>* resEdges) {
  const RouteEdgeInit<N, E, C>* curEdge = &settled.find(curE)->second;
  if (resNodes) resNodes->push_back(curEdge->e->getOtherNd(curEdge->n));
  while (resNodes || resEdges) {
    if (resNodes && curEdge->n) resNodes->push_back(curEdge->n);
    if (resEdges) resEdges->push_back(curEdge->e);
    if (!curEdge->parent) break;
    curEdge = &settled.find(curEdge->parent)->second;
  }
}
