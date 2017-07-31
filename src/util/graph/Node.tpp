// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>::Node(const N& pl) : _pl(pl) {
}

// _____________________________________________________________________________
template <typename N, typename E>
Node<N, E>::~Node() {
  for (auto e = _adjListOut.begin(); e != _adjListOut.end();) {
    Edge<N, E>* eP = *e;

    if (eP->getFrom() == this) {
      // careful with invalidating iterators
      e = _adjListOut.erase(e);
    } else {
      eP->getFrom()->removeEdge(eP);
      e++;
    }

    eP->getTo()->removeEdge(eP);

    delete eP;
  }

  for (auto e = _adjListIn.begin(); e != _adjListIn.end();) {
    Edge<N, E>* eP = *e;

    if (eP->getTo() == this) {
      // careful with invalidating iterators
      e = _adjListIn.erase(e);
    } else {
      eP->getTo()->removeEdge(eP);
      e++;
    }

    eP->getFrom()->removeEdge(eP);

    delete eP;
  }
}

// _____________________________________________________________________________
template <typename N, typename E>
void Node<N, E>::addEdge(Edge<N, E>* e) {
  if (e->getFrom() == this) _adjListOut.insert(e);
  if (e->getTo() == this) _adjListIn.insert(e);
}

// _____________________________________________________________________________
template <typename N, typename E>
void Node<N, E>::removeEdge(Edge<N, E>* e) {
  if (e->getFrom() == this) _adjListOut.erase(e);
  if (e->getTo() == this) _adjListIn.erase(e);
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::set<Edge<N, E>*>& Node<N, E>::getAdjListOut() const {
  return _adjListOut;
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::set<Edge<N, E>*>& Node<N, E>::getAdjListIn() const {
  return _adjListIn;
}

// _____________________________________________________________________________
template <typename N, typename E>
std::set<Edge<N, E>*> Node<N, E>::getAdjList() const {
  std::set<Edge<N, E>*> edges;
  edges.insert(getAdjListIn().begin(), getAdjListIn().end());
  edges.insert(getAdjListOut().begin(), getAdjListOut().end());

  return edges;
}

// _____________________________________________________________________________
template <typename N, typename E>
N& Node<N, E>::pl() {
  return _pl;
}

// _____________________________________________________________________________
template <typename N, typename E>
const N& Node<N, E>::pl() const {
  return _pl;
}