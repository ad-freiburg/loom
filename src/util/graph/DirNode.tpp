// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename N, typename E>
DirNode<N, E>::DirNode() : _pl() {}

// _____________________________________________________________________________
template <typename N, typename E>
DirNode<N, E>::DirNode(const N& pl) : _pl(pl) {}

// _____________________________________________________________________________
template <typename N, typename E>
DirNode<N, E>::~DirNode() {
  // delete self edges
  for (auto e = _adjListOut.begin(); e != _adjListOut.end();) {
    Edge<N, E>* eP = *e;
    if (eP->getTo() == this) {
      _adjListIn.erase(std::find(_adjListIn.begin(), _adjListIn.end(), eP));
      e = _adjListOut.erase(e);
      delete eP;
    } else {
      e++;
    }
  }

  for (auto e = _adjListOut.begin(); e != _adjListOut.end(); e++) {
    Edge<N, E>* eP = *e;

    if (eP->getTo() != this) {
      eP->getTo()->removeEdge(eP);
      delete eP;
    }
  }

  for (auto e = _adjListIn.begin(); e != _adjListIn.end(); e++) {
    Edge<N, E>* eP = *e;

    if (eP->getFrom() != this) {
      eP->getFrom()->removeEdge(eP);
      delete eP;
    }
  }
}

// _____________________________________________________________________________
template <typename N, typename E>
void DirNode<N, E>::addEdge(Edge<N, E>* e) {
  if (e->getFrom() == this && !adjOutContains(e)) {
    _adjListOut.reserve(_adjListOut.size() + 1);
    _adjListOut.push_back(e);
  }
  if (e->getTo() == this && !adjInContains(e)) {
    _adjListIn.reserve(_adjListIn.size() + 1);
    _adjListIn.push_back(e);
  }
}

// _____________________________________________________________________________
template <typename N, typename E>
void DirNode<N, E>::removeEdge(Edge<N, E>* e) {
  if (e->getFrom() == this) {
    auto p = std::find(_adjListOut.begin(), _adjListOut.end(), e);
    if (p != _adjListOut.end()) _adjListOut.erase(p);
  }
  if (e->getTo() == this) {
    auto p = std::find(_adjListIn.begin(), _adjListIn.end(), e);
    if (p != _adjListIn.end()) _adjListIn.erase(p);
  }
}
//
// _____________________________________________________________________________
template <typename N, typename E>
bool DirNode<N, E>::hasEdgeIn(const Edge<N, E>* e) const {
  return e->getTo() == this;
}

// _____________________________________________________________________________
template <typename N, typename E>
bool DirNode<N, E>::hasEdgeOut(const Edge<N, E>* e) const {
  return e->getFrom() == this;
}

// _____________________________________________________________________________
template <typename N, typename E>
bool DirNode<N, E>::hasEdge(const Edge<N, E>* e) const {
  return hasEdgeOut(e) || hasEdgeIn(e);
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::vector<Edge<N, E>*>& DirNode<N, E>::getAdjList() const {
  return _adjListOut;
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::vector<Edge<N, E>*>& DirNode<N, E>::getAdjListOut() const {
  return _adjListOut;
}

// _____________________________________________________________________________
template <typename N, typename E>
const std::vector<Edge<N, E>*>& DirNode<N, E>::getAdjListIn() const {
  return _adjListIn;
}

// _____________________________________________________________________________
template <typename N, typename E>
size_t DirNode<N, E>::getDeg() const {
  return _adjListOut.size();
}

// _____________________________________________________________________________
template <typename N, typename E>
size_t DirNode<N, E>::getInDeg() const {
  return _adjListIn.size();
}

// _____________________________________________________________________________
template <typename N, typename E>
size_t DirNode<N, E>::getOutDeg() const {
  return _adjListOut.size();
}

// _____________________________________________________________________________
template <typename N, typename E>
N& DirNode<N, E>::pl() {
  return _pl;
}

// _____________________________________________________________________________
template <typename N, typename E>
const N& DirNode<N, E>::pl() const {
  return _pl;
}

// _____________________________________________________________________________
template <typename N, typename E>
bool DirNode<N, E>::adjInContains(const Edge<N, E>* e) const {
  for (size_t i = 0; i < _adjListIn.size(); i++)
    if (_adjListIn[i] == e) return true;
  return false;
}

// _____________________________________________________________________________
template <typename N, typename E>
bool DirNode<N, E>::adjOutContains(const Edge<N, E>* e) const {
  for (size_t i = 0; i < _adjListOut.size(); i++)
    if (_adjListOut[i] == e) return true;
  return false;
}
