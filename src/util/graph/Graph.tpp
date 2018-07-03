// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
Graph<N, E>::~Graph() {
  for (auto n : _nodes) delete n;
}

// _____________________________________________________________________________
template <typename N, typename E>
Edge<N, E>* Graph<N, E>::addEdg(Node<N, E>* from, Node<N, E>* to) {
  return addEdg(from, to, E());
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::set<Node<N, E>*>& Graph<N, E>::getNds() const {
  return _nodes;
}

// _____________________________________________________________________________
template <typename N, typename E>
std::set<Node<N, E>*>* Graph<N, E>::getNds() {
  return &_nodes;
}

// _____________________________________________________________________________
template <typename N, typename E>
typename std::set<Node<N, E>*>::iterator Graph<N, E>::delNd(
    Node<N, E>* n) {
  return delNd(_nodes.find(n));
}

// _____________________________________________________________________________
template <typename N, typename E>
typename std::set<Node<N, E>*>::iterator Graph<N, E>::delNd(
    typename std::set<Node<N, E>*>::iterator i) {
  delete *i;
  return _nodes.erase(i);
}

// _____________________________________________________________________________
template <typename N, typename E>
void Graph<N, E>::delEdg(Node<N, E>* from, Node<N, E>* to) {
  Edge<N, E>* toDel = getEdg(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdg(from, to));

  delete toDel;
}

// _____________________________________________________________________________
template <typename N, typename E>
Edge<N, E>* Graph<N, E>::getEdg(Node<N, E>* from, Node<N, E>* to) {
  for (auto e : from->getAdjList()) {
    if (e->getOtherNd(from) == to) return e;
  }

  return 0;
}

// _____________________________________________________________________________
template <typename N, typename E>
const Edge<N, E>* Graph<N, E>::getEdg(Node<N, E>* from, Node<N, E>* to) const {
  for (auto e : from->getAdjList()) {
    if (e->getOtherNd(from) == to) return e;
  }

  return 0;
}
