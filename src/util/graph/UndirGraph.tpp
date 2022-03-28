// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
UndirGraph<N, E>::UndirGraph() {}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* UndirGraph<N, E>::addNd(const N& pl) {
  return addNd(new UndirNode<N, E>(pl));
}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* UndirGraph<N, E>::addNd() {
  return addNd(new UndirNode<N, E>());
}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* UndirGraph<N, E>::addNd(UndirNode<N, E>* n) {
  auto ins = Graph<N, E>::_nodes.insert(n);
  return *ins.first;
}

// _____________________________________________________________________________
template <typename N, typename E>
Edge<N, E>* UndirGraph<N, E>::addEdg(Node<N, E>* from, Node<N, E>* to,
                                     const E& p) {
  Edge<N, E>* e = Graph<N, E>::getEdg(from, to);
  if (!e) {
    e = new Edge<N, E>(from, to, p);
    from->addEdge(e);
    to->addEdge(e);
  }
  return e;
}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* UndirGraph<N, E>::mergeNds(Node<N, E>* a, Node<N, E>* b) {
  for (auto e : a->getAdjListOut()) {
    if (e->getFrom() != a) continue;
    if (e->getTo() != b) {
      addEdg(b, e->getTo(), e->pl());
    }
  }

  for (auto e : a->getAdjListIn()) {
    if (e->getTo() != a) continue;
    if (e->getFrom() != b) {
      addEdg(e->getFrom(), b, e->pl());
    }
  }

  UndirGraph<N, E>::delNd(a);

  return b;
}
