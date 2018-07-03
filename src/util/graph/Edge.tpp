// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
Edge<N, E>::Edge(Node<N, E>* from, Node<N, E>* to, const E& pl)
 : _from(from), _to(to), _pl(pl) {

}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* Edge<N, E>::getFrom() const {
  return _from;
}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* Edge<N, E>::getTo() const {
  return _to;
}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>* Edge<N, E>::getOtherNd(const Node<N, E>* notNode) const {
  if (_to == notNode) return _from;
  return _to;
}

// _____________________________________________________________________________
template <typename N, typename E>
E& Edge<N, E>::pl() {
  return _pl;
}

// _____________________________________________________________________________
template <typename N, typename E>
const E& Edge<N, E>::pl() const {
  return _pl;
}
