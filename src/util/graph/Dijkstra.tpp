// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
bool Dijkstra::shortestPath(Node<N, E>* from,
                            const std::unordered_map<Node<N, E>*, bool>& to,
                            std::list<Edge<N, E>*>* res, Node<N, E>** target) {
  std::unordered_map<Node<N, E>*, RouteNode<N, E> > settled;
  std::priority_queue<RouteNode<N, E> > pq;
  bool found = false;

  RouteNode<N, E> start(from, 0, 0, 0);
  pq.push(start);
  RouteNode<N, E> cur = start;

  while (!pq.empty()) {
    cur = pq.top();
    pq.pop();

    if (settled.find(cur.n) != settled.end()) continue;
    settled[cur.n] = cur;

    if (to.find(cur.n) != to.end()) {
      found = true;
      break;
    }

    // relaxation
    for (auto edge : cur.n->getAdjListOut()) {
      double newC = cur.d + edge->pl().cost();
      pq.push(RouteNode<N, E>(edge->getTo(), cur.n, newC, &(*edge)));
    }
  }

  if (!found) return false;

  // traverse the parents backwards beginning at current target node
  Node<N, E>* curN = cur.n;
  *target = cur.n;

  while (true) {
    const RouteNode<N, E>& curNode = settled[curN];
    if (curN == from) break;
    res->push_front(curNode.e);
    curN = curNode.parent;
  }

  return true;
}

// _____________________________________________________________________________
template <typename N, typename E>
bool Dijkstra::shortestPath(Node<N, E>* from, Node<N, E>* to,
                            std::list<Edge<N, E>*>* res) {
  Node<N, E>* target;
  std::unordered_map<Node<N, E>*, bool> tos;
  tos[to] = true;
  return shortestPath(from, tos, res, &target);
}
