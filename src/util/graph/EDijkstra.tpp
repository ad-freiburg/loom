// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C EDijkstra::shortestPathImpl(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                              const ShortestPath::CostFunc<N, E, C>& costFunc,
                              const ShortestPath::HeurFunc<N, E, C>& heurFunc,
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
                              const ShortestPath::CostFunc<N, E, C>& costFunc,
                              const ShortestPath::HeurFunc<N, E, C>& heurFunc,
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
C EDijkstra::shortestPathImpl(const std::set<Edge<N, E>*> from,
                              const std::set<Edge<N, E>*>& to,
                              const ShortestPath::CostFunc<N, E, C>& costFunc,
                              const ShortestPath::HeurFunc<N, E, C>& heurFunc,
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
    pq.emplace(e, (Edge<N, E>*)0, (Node<N, E>*)0, c, c + h);
  }

  RouteEdge<N, E, C> cur;

  while (!pq.empty()) {
    EDijkstra::ITERS++;

    if (settled.find(pq.top().e) != settled.end()) {
      pq.pop();
      continue;
    }

    cur = pq.top();
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
    const std::set<Edge<N, E>*>& from, const ShortestPath::CostFunc<N, E, C>& costFunc,
    bool rev) {
  std::unordered_map<Edge<N, E>*, C> costs;

  Settled<N, E, C> settled;
  PQ<N, E, C> pq;

  std::set<Edge<N, E>*> to;

  for (auto e : from) {
    pq.emplace(e, (Edge<N, E>*)0, (Node<N, E>*)0, costFunc(0, 0, e), C());
  }

  RouteEdge<N, E, C> cur;

  while (!pq.empty()) {
    EDijkstra::ITERS++;

    if (settled.find(pq.top().e) != settled.end()) {
      pq.pop();
      continue;
    }

    cur = pq.top();
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
    const ShortestPath::CostFunc<N, E, C>& costFunc,
    const ShortestPath::HeurFunc<N, E, C>& heurFunc,
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
  pq.emplace(from, (Edge<N, E>*)0, (Node<N, E>*)0, c, c + h);

  RouteEdge<N, E, C> cur;

  while (!pq.empty()) {
    EDijkstra::ITERS++;

    if (settled.find(pq.top().e) != settled.end()) {
      pq.pop();
      continue;
    }

    cur = pq.top();
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
void EDijkstra::relaxInv(RouteEdge<N, E, C>& cur,
                         const ShortestPath::CostFunc<N, E, C>& costFunc,
                         PQ<N, E, C>& pq) {

  // handling undirected graph makes no sense here

  for (const auto edge : cur.e->getFrom()->getAdjListIn()) {
    if (edge == cur.e) continue;
    C newC = costFunc(edge, cur.e->getFrom(), cur.e);
    newC = cur.d + newC;
    if (costFunc.inf() <= newC) continue;

    pq.emplace(edge, cur.e, cur.e->getFrom(), newC, C());
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::relax(RouteEdge<N, E, C>& cur, const std::set<Edge<N, E>*>& to,
                      const ShortestPath::CostFunc<N, E, C>& costFunc,
                      const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                      PQ<N, E, C>& pq) {
  if (cur.e->getFrom()->hasEdgeIn(cur.e)) {
    // for undirected graphs
    for (const auto edge : cur.e->getFrom()->getAdjListOut()) {
      if (edge == cur.e) continue;
      C newC = costFunc(cur.e, cur.e->getFrom(), edge);
      newC = cur.d + newC;
      if (costFunc.inf() <= newC) continue;

      const C& h = heurFunc(edge, to);
      const C& newH = newC + h;

      pq.emplace(edge, cur.e, cur.e->getFrom(), newC, newH);
    }
  }

  for (const auto edge : cur.e->getTo()->getAdjListOut()) {
    if (edge == cur.e) continue;
    C newC = costFunc(cur.e, cur.e->getTo(), edge);
    newC = cur.d + newC;
    if (costFunc.inf() <= newC) continue;

    const C& h = heurFunc(edge, to);
    const C& newH = newC + h;

    pq.emplace(edge, cur.e, cur.e->getTo(), newC, newH);
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void EDijkstra::buildPath(Edge<N, E>* curE, const Settled<N, E, C>& settled,
                          NList<N, E>* resNodes, EList<N, E>* resEdges) {
  const RouteEdge<N, E, C>* curEdge = &settled.find(curE)->second;
  if (resNodes) resNodes->push_back(curEdge->e->getOtherNd(curEdge->n));
  while (true) {
    if (resNodes && curEdge->n) resNodes->push_back(curEdge->n);
    if (resEdges) resEdges->push_back(curEdge->e);
    if (!curEdge->parent) break;
    curEdge = &settled.find(curEdge->parent)->second;
  }
}
