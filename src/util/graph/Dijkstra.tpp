// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
int Dijkstra::shortestPath(Node<N, E>* from,
                            const std::unordered_set<Node<N, E>*>& to,
                            std::list<Edge<N, E>*>* res, Node<N, E>** target) {
  std::unordered_map<Node<N, E>*, RouteNode<N, E> > settled;
  std::priority_queue<RouteNode<N, E> > pq;
  bool found = false;
  int iter = 0;

  RouteNode<N, E> start(from, 0, 0, 0);
  pq.push(start);
  RouteNode<N, E> cur = start;

  while (!pq.empty()) {
    iter++;
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
      if (edge->pl().cost() == std::numeric_limits<double>::infinity()) continue;
      double newC = cur.d + edge->pl().cost();
      pq.push(RouteNode<N, E>(edge->getTo(), cur.n, newC, &(*edge)));
    }
  }

  if (!found) return 0;

  // traverse the parents backwards beginning at current target node
  Node<N, E>* curN = cur.n;
  *target = cur.n;

  while (true) {
    const RouteNode<N, E>& curNode = settled[curN];
    if (curN == from) break;
    res->push_front(curNode.e);
    curN = curNode.parent;
  }

  return iter;
}

// _____________________________________________________________________________
template <typename N, typename E, typename H>
int Dijkstra::shortestPathAStar(Node<N, E>* from,
                            const std::unordered_set<Node<N, E>*>& to,
                            H heurFunc,
                            std::list<Edge<N, E>*>* res, Node<N, E>** target,
                            bool check) {
  std::unordered_map<Node<N, E>*, RouteNodeAStar<N, E> > settled;
  std::priority_queue<RouteNodeAStar<N, E> > pq;
  bool found = false;
  int iter = 0;

  RouteNodeAStar<N, E> start(from, 0, 0, 0, 0);
  pq.push(start);
  RouteNodeAStar<N, E> cur = start;

  while (!pq.empty()) {
    iter++;
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
      if (edge->pl().cost() == std::numeric_limits<double>::infinity()) continue;
      double newC = cur.d + edge->pl().cost();

      double hCost = heurFunc(edge->getTo(), to);
      double newH = newC + hCost;

      if(check) edge->pl().setVisited(iter);

      pq.push(RouteNodeAStar<N, E>(edge->getTo(), cur.n, newH, newC, &(*edge)));
    }
  }

  if (!found) return 0;

  // traverse the parents backwards beginning at current target node
  Node<N, E>* curN = cur.n;
  *target = cur.n;

  while (true) {
    const RouteNodeAStar<N, E>& curNode = settled[curN];
    if (curN == from) break;
    res->push_front(curNode.e);
    curN = curNode.parent;
  }

  return iter;
}

// _____________________________________________________________________________
template <typename N, typename E>
int Dijkstra::shortestPath(Node<N, E>* from, Node<N, E>* to,
                            std::list<Edge<N, E>*>* res) {
  Node<N, E>* target;
  std::unordered_set<Node<N, E>*, bool> tos;
  tos.insert(to);
  return shortestPath(from, tos, res, &target);
}

// _____________________________________________________________________________
template <typename N, typename E, typename H>
int Dijkstra::shortestPathAStar(Node<N, E>* from, Node<N, E>* to,
                            H heur, std::list<Edge<N, E>*>* res) {
  Node<N, E>* target;
  std::unordered_set<Node<N, E>*, bool> tos;
  tos.insert(to);
  return shortestPathAStar(from, tos, res, heur, &target, false);
}
