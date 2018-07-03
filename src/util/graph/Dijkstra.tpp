// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C Dijkstra::shortestPathImpl(Node<N, E>* from, const std::set<Node<N, E>*>& to,
                             const ShortestPath::CostFunc<N, E, C>& costFunc,
                             const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                             EList<N, E>* resEdges, NList<N, E>* resNodes) {
  if (from->getOutDeg() == 0) return costFunc.inf();

  Settled<N, E, C> settled;
  PQ<N, E, C> pq;
  bool found = false;

  pq.emplace(from);
  RouteNode<N, E, C> cur;

  while (!pq.empty()) {
    Dijkstra::ITERS++;

    if (settled.find(pq.top().n) != settled.end()) {
      pq.pop();
      continue;
    }

    cur = pq.top();
    pq.pop();

    settled[cur.n] = cur;

    if (to.find(cur.n) != to.end()) {
      found = true;
      break;
    }

    relax(cur, to, costFunc, heurFunc, pq);
  }

  if (!found) return costFunc.inf();

  buildPath(cur.n, settled, resNodes, resEdges);

  return cur.d;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C Dijkstra::shortestPathImpl(const std::set<Node<N, E>*> from,
                             const std::set<Node<N, E>*>& to,
                             const ShortestPath::CostFunc<N, E, C>& costFunc,
                             const ShortestPath::HeurFunc<N, E, C>& heurFunc,
                             EList<N, E>* resEdges, NList<N, E>* resNodes) {
  Settled<N, E, C> settled;
  PQ<N, E, C> pq;
  bool found = false;

  // put all nodes in from onto PQ
  for (auto n : from) pq.emplace(n);
  RouteNode<N, E, C> cur;

  while (!pq.empty()) {
    Dijkstra::ITERS++;

    if (settled.find(pq.top().n) != settled.end()) {
      pq.pop();
      continue;
    }

    cur = pq.top();
    pq.pop();

    settled[cur.n] = cur;

    if (to.find(cur.n) != to.end()) {
      found = true;
      break;
    }

    relax(cur, to, costFunc, heurFunc, pq);
  }

  if (!found) return costFunc.inf();

  buildPath(cur.n, settled, resNodes, resEdges);

  return cur.d;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
std::unordered_map<Node<N, E>*, C> Dijkstra::shortestPathImpl(
    Node<N, E>* from, const std::set<Node<N, E>*>& to,
    const ShortestPath::CostFunc<N, E, C>& costFunc,
    const ShortestPath::HeurFunc<N, E, C>& heurFunc,
    std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges,
    std::unordered_map<Node<N, E>*, NList<N, E>*> resNodes) {
  std::unordered_map<Node<N, E>*, C> costs;
  if (to.size() == 0) return costs;
  // init costs with inf
  for (auto n : to) costs[n] = costFunc.inf();

  if (from->getOutDeg() == 0) return costs;

  Settled<N, E, C> settled;
  PQ<N, E, C> pq;

  size_t found = 0;

  pq.emplace(from);
  RouteNode<N, E, C> cur;

  while (!pq.empty()) {
    Dijkstra::ITERS++;

    if (settled.find(pq.top().n) != settled.end()) {
      pq.pop();
      continue;
    }

    cur = pq.top();
    pq.pop();

    settled[cur.n] = cur;

    if (to.find(cur.n) != to.end()) {
      found++;
    }

    if (found == to.size()) break;

    relax(cur, to, costFunc, heurFunc, pq);
  }

  for (auto nto : to) {
    if (!settled.count(nto)) continue;
    Node<N, E>* curN = nto;
    costs[nto] = settled[curN].d;

    buildPath(nto, settled, resNodes[nto], resEdges[nto]);
  }

  return costs;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void Dijkstra::relax(RouteNode<N, E, C>& cur, const std::set<Node<N, E>*>& to,
                     const ShortestPath::CostFunc<N, E, C>& costFunc,
                     const ShortestPath::HeurFunc<N, E, C>& heurFunc, PQ<N, E, C>& pq) {
  for (auto edge : cur.n->getAdjListOut()) {
    C newC = costFunc(cur.n, edge, edge->getOtherNd(cur.n));
    newC = cur.d + newC;
    if (costFunc.inf() <= newC) continue;

    const C& newH = newC + heurFunc(edge->getOtherNd(cur.n), to);

    pq.emplace(edge->getOtherNd(cur.n), cur.n, newC, newH, &(*edge));
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void Dijkstra::buildPath(Node<N, E>* curN,
                         Settled<N, E, C>& settled, NList<N, E>* resNodes,
                         EList<N, E>* resEdges) {
  while (true) {
    const RouteNode<N, E, C>& curNode = settled[curN];
    if (resNodes) resNodes->push_back(curNode.n);
    if (!curNode.parent) break;
    if (resEdges) resEdges->push_back(curNode.e);
    curN = curNode.parent;
  }
}
