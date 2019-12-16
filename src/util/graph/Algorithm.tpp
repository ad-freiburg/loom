// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>
//
// _____________________________________________________________________________
template <typename N, typename E>
std::vector<std::set<Node<N, E>*>> Algorithm::connectedComponents(
    const UndirGraph<N, E>& g) {
  return connectedComponents(g, EdgeCheckFunc<N, E>());
}

// _____________________________________________________________________________
template <typename N, typename E>
std::vector<std::set<Node<N, E>*>> Algorithm::connectedComponents(
    const UndirGraph<N, E>& g, const EdgeCheckFunc<N, E>& checkFunc) {
  std::vector<std::set<Node<N, E>*>> ret;
  std::set<Node<N, E>*> visited;

  for (auto* n : g.getNds()) {
    if (!visited.count(n)) {
      ret.resize(ret.size() + 1);
      std::stack<Node<N, E>*> q;
      q.push(n);
      while (!q.empty()) {
        Node<N, E>* cur = q.top();
        q.pop();

        ret.back().insert(cur);
        visited.insert(cur);

        for (auto* e : cur->getAdjList()) {
          if (!checkFunc(cur, e)) continue;
          if (!visited.count(e->getOtherNd(cur))) q.push(e->getOtherNd(cur));
        }
      }
    }
  }

  return ret;
}
