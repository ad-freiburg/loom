// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C BiDijkstra::shortestPathImpl(Node<N, E>* from,
                               const std::set<Node<N, E>*>& to,
                               const util::graph::CostFunc<N, E, C>& costFunc,
                               const util::graph::HeurFunc<N, E, C>& heurFunc,
                               EList<N, E>* resEdges, NList<N, E>* resNodes) {
  if (from->getOutDeg() == 0) return costFunc.inf();

  std::set<Node<N, E>*> froms;
  froms.insert(from);

  return shortestPathImpl(froms, to, costFunc, heurFunc, resEdges, resNodes);
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C BiDijkstra::shortestPathImpl(const std::set<Node<N, E>*> from,
                               const std::set<Node<N, E>*>& to,
                               const util::graph::CostFunc<N, E, C>& costFunc,
                               const util::graph::HeurFunc<N, E, C>& heurFunc,
                               EList<N, E>* resEdges, NList<N, E>* resNodes) {
  Settled<N, E, C> settledFwd, settledBwd;
  PQ<N, E, C> pqFwd, pqBwd;
  bool found = false;

  // starter for forward search
  for (auto n : from) pqFwd.emplace(n);

  auto l = costFunc.inf();

  // starter for backward search
  for (auto n : to) pqBwd.emplace(n);

  RouteNode<N, E, C> cur;

  while (!pqFwd.empty() && !pqBwd.empty()) {
    if (costFunc.inf() <= pqFwd.top().h && costFunc.inf() <= pqBwd.top().h)
      return costFunc.inf();

    if (pqFwd.top() < pqBwd.top()) {
      auto se = settledBwd.find(pqBwd.top().n);
      if (se != settledBwd.end()) {
        // to allow non-consistent heuristics
        if (se->second.d <= pqBwd.top().d) {
          pqBwd.pop();
          continue;
        }
      }
    } else {
      auto se = settledFwd.find(pqFwd.top().n);
      if (se != settledFwd.end()) {
        // to allow non-consistent heuristics
        if (se->second.d <= pqFwd.top().d) {
          pqFwd.pop();
          continue;
        }
      }
    }

    BiDijkstra::ITERS++;

    if (pqFwd.top() < pqBwd.top()) {
      cur = pqBwd.top();
      pqBwd.pop();
      settledBwd[cur.n] = cur;
      if (settledFwd.find(cur.n) != settledFwd.end()) {
        auto newL = cur.d + settledFwd.find(cur.n)->second.d;

        if (!(newL > l)) {
          l = newL;
          found = true;
          break;
        }
      }
      C bestCost = relaxBwd(from, cur, costFunc, heurFunc, pqBwd, settledFwd);
      if (bestCost < l) l = bestCost;
    } else {
      cur = pqFwd.top();
      pqFwd.pop();
      settledFwd[cur.n] = cur;
      if (settledBwd.find(cur.n) != settledBwd.end()) {
        auto newL = cur.d + settledBwd.find(cur.n)->second.d;

        if (!(newL > l)) {
          l = newL;
          found = true;
          break;
        }
      }
      C bestCost = relaxFwd(cur, to, costFunc, heurFunc, pqFwd, settledBwd);
      if (bestCost < l) l = bestCost;
    }
  }

  if (!found) return costFunc.inf();

  buildPath(cur.n, settledFwd, settledBwd, resNodes, resEdges);

  return l;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
std::unordered_map<Node<N, E>*, C> BiDijkstra::shortestPathImpl(
    Node<N, E>* from, const std::set<Node<N, E>*>& to,
    const util::graph::CostFunc<N, E, C>& costFunc,
    const util::graph::HeurFunc<N, E, C>& heurFunc,
    std::unordered_map<Node<N, E>*, EList<N, E>*> resEdges,
    std::unordered_map<Node<N, E>*, NList<N, E>*> resNodes) {
  UNUSED(from);
  UNUSED(to);
  UNUSED(costFunc);
  UNUSED(heurFunc);
  UNUSED(resEdges);
  UNUSED(resNodes);
  assert(false);
  // std::unordered_map<Node<N, E>*, C> costs;
  // if (to.size() == 0) return costs;
  // // init costs with inf
  // for (auto n : to) costs[n] = costFunc.inf();

  // if (from->getOutDeg() == 0) return costs;
  // Settled<N, E, C> settled;
  // PQ<N, E, C> pq;

  // size_t found = 0;

  // pq.emplace(from);
  // RouteNode<N, E, C> cur;

  // while (!pq.empty()) {
  // if (costFunc.inf() <= pq.top().h) return costs;
  // if (settled.find(pq.top().n) != settled.end()) {
  // pq.pop();
  // continue;
  // }
  // BiDijkstra::ITERS++;

  // cur = pq.top();
  // pq.pop();

  // settled[cur.n] = cur;

  // if (to.find(cur.n) != to.end()) {
  // found++;
  // }

  // if (found == to.size()) break;

  // relax(cur, to, costFunc, heurFunc, pq);
  // }

  // for (auto nto : to) {
  // if (!settled.count(nto)) continue;
  // Node<N, E>* curN = nto;
  // costs[nto] = settled[curN].d;

  // buildPath(nto, settled, resNodes[nto], resEdges[nto]);
  // }

  // return costs;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void BiDijkstra::relax(RouteNode<N, E, C>& cur, const std::set<Node<N, E>*>& to,
                       const util::graph::CostFunc<N, E, C>& costFunc,
                       const util::graph::HeurFunc<N, E, C>& heurFunc,
                       PQ<N, E, C>& pq) {
  for (auto edge : cur.n->getAdjListOut()) {
    C newC = costFunc(cur.n, edge, edge->getOtherNd(cur.n));
    newC = cur.d + newC;
    if (costFunc.inf() <= newC) continue;

    // addition done here to avoid it in the PQ
    const C& newH = newC + heurFunc(edge->getOtherNd(cur.n), to);

    pq.emplace(edge->getOtherNd(cur.n), cur.n, newC, newH);
  }
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C BiDijkstra::relaxFwd(RouteNode<N, E, C>& cur, const std::set<Node<N, E>*>& to,
                       const util::graph::CostFunc<N, E, C>& costFunc,
                       const util::graph::HeurFunc<N, E, C>& heurFunc,
                       PQ<N, E, C>& pq, const Settled<N, E, C>& settledBwd) {
  UNUSED(to);
  UNUSED(heurFunc);
  C ret = costFunc.inf();
  for (auto edge : cur.n->getAdjListOut()) {
    C newC = costFunc(cur.n, edge, edge->getOtherNd(cur.n));
    newC = cur.d + newC;
    if (costFunc.inf() <= newC) continue;

    // addition done here to avoid it in the PQ
    // const C& newH = newC + heurFunc(froms, edge->getOtherNd(cur.n));

    // TODO:
    const C& newH = newC + 0;

    // update new best found cost
    if (settledBwd.find(edge->getOtherNd(cur.n)) != settledBwd.end()) {
      C bwdCost = settledBwd.find(edge->getOtherNd(cur.n))->second.d + newC;
      if (bwdCost < ret) ret = bwdCost;
    }

    pq.emplace(edge->getOtherNd(cur.n), cur.n, newC, newH);
  }

  return ret;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
C BiDijkstra::relaxBwd(const std::set<Node<N, E>*>& froms,
                       RouteNode<N, E, C>& cur,
                       const util::graph::CostFunc<N, E, C>& costFunc,
                       const util::graph::HeurFunc<N, E, C>& heurFunc,
                       PQ<N, E, C>& pq, const Settled<N, E, C>& settledFwd) {
  UNUSED(froms);
  UNUSED(heurFunc);
  C ret = costFunc.inf();
  for (auto edge : cur.n->getAdjListIn()) {
    C newC = costFunc(edge->getOtherNd(cur.n), edge, cur.n);
    newC = cur.d + newC;
    if (costFunc.inf() <= newC) continue;

    // addition done here to avoid it in the PQ
    // const C& newH = newC + heurFunc(froms, edge->getOtherNd(cur.n));

    // TODO:
    const C& newH = newC + 0;

    // update new best found cost
    if (settledFwd.find(edge->getOtherNd(cur.n)) != settledFwd.end()) {
      C fwdCost = settledFwd.find(edge->getOtherNd(cur.n))->second.d + newC;
      if (fwdCost < ret) ret = fwdCost;
    }

    pq.emplace(edge->getOtherNd(cur.n), cur.n, newC, newH);
  }

  return ret;
}

// _____________________________________________________________________________
template <typename N, typename E, typename C>
void BiDijkstra::buildPath(Node<N, E>* curN, Settled<N, E, C>& settledFwd,
                           Settled<N, E, C>& settledBwd, NList<N, E>* resNodes,
                           EList<N, E>* resEdges) {
  Node<N, E>* curNFwd = curN;
  Node<N, E>* curNBwd = curN;

  // the forward part
  while (resNodes || resEdges) {
    const RouteNode<N, E, C>& curNode = settledFwd[curNFwd];
    if (resNodes) resNodes->push_back(curNode.n);
    if (!curNode.parent) break;

    if (resEdges) {
      for (auto e : curNode.n->getAdjListIn()) {
        if (e->getOtherNd(curNode.n) == curNode.parent) resEdges->push_back(e);
      }
    }
    curNFwd = curNode.parent;
  }

  if (resNodes) std::reverse(resNodes->begin(), resNodes->end());
  if (resEdges) std::reverse(resEdges->begin(), resEdges->end());

  // the backward part
  while (resNodes || resEdges) {
    const RouteNode<N, E, C>& curNode = settledBwd[curNBwd];
    if (resNodes && curNode.n != curN) resNodes->push_back(curNode.n);
    if (!curNode.parent) break;

    if (resEdges) {
      for (auto e : curNode.n->getAdjListOut()) {
        if (e->getOtherNd(curNode.n) == curNode.parent) resEdges->push_back(e);
      }
    }
    curNBwd = curNode.parent;
  }

  if (resNodes) std::reverse(resNodes->begin(), resNodes->end());
  if (resEdges) std::reverse(resEdges->begin(), resEdges->end());
}
