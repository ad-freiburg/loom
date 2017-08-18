// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
Graph<N, E>::Graph() {
}

// _____________________________________________________________________________
template <typename N, typename E>
Graph<N, E>::~Graph() {
  for (auto n : _nodes) {
    delete n;
  }
}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* Graph<N, E>::addNode(Node<N, E>* n) {
  auto ins = _nodes.insert(n);
  return *ins.first;
}

// _____________________________________________________________________________
template <typename N, typename E>
Edge<N, E>* Graph<N, E>::addEdge(Node<N, E>* from, Node<N, E>* to, const E& p) {
  if (from == to) return 0;
  Edge<N, E>* e = getEdge(from, to);
  if (!e) {
    e = new Edge<N, E>(from, to, p);
    from->addEdge(e);
    to->addEdge(e);
  }
  return e;
}

// _____________________________________________________________________________
template <typename N, typename E>
Edge<N, E>* Graph<N, E>::getEdge(Node<N, E>* from, Node<N, E>* to) {
  for (auto e : from->getAdjListOut()) {
    if (e->getTo() == to) return e;
  }

  // also search in the opposite direction, we are handling an undirected
  // graph here
  for (auto e : from->getAdjListIn()) {
    if (e->getFrom() == to) return e;
  }

  return 0;
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::set<Node<N, E>*>& Graph<N, E>::getNodes() const { return _nodes; }

// _____________________________________________________________________________
template <typename N, typename E>
std::set<Node<N, E>*>* Graph<N, E>::getNodes() { return &_nodes; }

// _____________________________________________________________________________
template <typename N, typename E>
void Graph<N, E>::deleteEdge(Node<N, E>* from, Node<N, E>* to) {
  Edge<N, E>* toDel = getEdge(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdge(from, to));

  delete toDel;
}

// _____________________________________________________________________________
template <typename N, typename E>
void Graph<N, E>::deleteNode(Node<N, E>* n) {\
  _nodes.erase(n); 
  delete n;
}
