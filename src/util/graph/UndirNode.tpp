// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
UndirNode<N, E>::UndirNode() : _pl() {
}

// _____________________________________________________________________________
template <typename N, typename E>
UndirNode<N, E>::UndirNode(const N& pl) : _pl(pl) {
}

// _____________________________________________________________________________
template <typename N, typename E>
UndirNode<N, E>::~UndirNode() {
  // delete self edges
  for (auto e = _adjList.begin(); e != _adjList.end();) {
    Edge<N, E>* eP = *e;
    if (eP->getTo() == this && eP->getFrom() == this) {
      e = _adjList.erase(e);
      delete eP;
    } else {
      e++;
    }
  }

  for (auto e = _adjList.begin(); e != _adjList.end(); e++) {
    Edge<N, E>* eP = *e;

    if (eP->getTo() != this) {
      eP->getTo()->removeEdge(eP);
    }

    if (eP->getFrom() != this) {
      eP->getFrom()->removeEdge(eP);
    }

    delete eP;
  }
}

// _____________________________________________________________________________
template <typename N, typename E>
bool UndirNode<N, E>::hasEdgeIn(const Edge<N, E>* e) const {
  return hasEdge(e);
}

// _____________________________________________________________________________
template <typename N, typename E>
bool UndirNode<N, E>::hasEdgeOut(const Edge<N, E>* e) const {
  return hasEdge(e);
}

// _____________________________________________________________________________
template <typename N, typename E>
bool UndirNode<N, E>::hasEdge(const Edge<N, E>* e) const {
  return e->getFrom() == this || e->getTo() == this;
}

// _____________________________________________________________________________
template <typename N, typename E>
bool UndirNode<N, E>::adjContains(const Edge<N, E>* e) const {
  for (size_t i = 0; i < _adjList.size(); i++) if (_adjList[i] == e) return true;
  return false;
}

// _____________________________________________________________________________
template <typename N, typename E>
void UndirNode<N, E>::addEdge(Edge<N, E>* e) {
  if (adjContains(e)) return;
  _adjList.reserve(_adjList.size() + 1);
  _adjList.push_back(e);
}

// _____________________________________________________________________________
template <typename N, typename E>
void UndirNode<N, E>::removeEdge(Edge<N, E>* e) {
  auto p = std::find(_adjList.begin(), _adjList.end(), e);
  if (p != _adjList.end()) _adjList.erase(p);
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::vector<Edge<N, E>*>& UndirNode<N, E>::getAdjList() const {
  return _adjList;
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::vector<Edge<N, E>*>& UndirNode<N, E>::getAdjListOut() const {
  return _adjList;
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::vector<Edge<N, E>*>& UndirNode<N, E>::getAdjListIn() const {
  return _adjList;
}

// _____________________________________________________________________________
template <typename N, typename E>
size_t UndirNode<N, E>::getDeg() const {
  return _adjList.size();
}

// _____________________________________________________________________________
template <typename N, typename E>
size_t UndirNode<N, E>::getInDeg() const {
  return getDeg();
}

// _____________________________________________________________________________
template <typename N, typename E>
size_t UndirNode<N, E>::getOutDeg() const {
  return getDeg();
}

// _____________________________________________________________________________
template <typename N, typename E>
N& UndirNode<N, E>::pl() {
  return _pl;
}

// _____________________________________________________________________________
template <typename N, typename E>
const N& UndirNode<N, E>::pl() const {
  return _pl;
}
