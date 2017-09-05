// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
bool Dijkstra::shortestPath(Node<N, E>* from, Node<N, E>* to,
                                   std::list<Edge<N, E>*>* res) {
  std::unordered_map<Node<N, E>*, bool> settled;
  std::unordered_map<Node<N, E>*, RouteNode<N, E> > finishedNodes;
  std::priority_queue<RouteNode<N, E> > pq;
  bool found = false;

  RouteNode<N, E> start(from, 0, 0, 0);
  pq.push(start);
  RouteNode<N, E> cur = start;

  while (!pq.empty()) {
    cur = pq.top();
    pq.pop();

    if (settled[cur.n]) continue;

    settled[cur.n] = true;
    finishedNodes[cur.n] = cur;

    if (cur.n == to) {
      found = true;
      break;
    }

    // relaxation
    for (auto edge : cur.n->getAdjListOut()) {
      // TODO: real costs!
      double newC = cur.d + 1;

      pq.push(RouteNode<N, E>(edge->getTo(), cur.n, newC, &(*edge)));
    }
  }

  if (!found) return false;

  // traverse the parents backwards beginning at current target node
  Node<N, E>* curN = cur.n;

  while (true) {
    const RouteNode<N, E>& curNode = finishedNodes[curN];
    if (curN == from) break;
    res->push_front(curNode.e);
    curN = curNode.parent;
  }

  // found shortest path
  return true;
}
