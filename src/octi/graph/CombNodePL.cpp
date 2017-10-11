// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/CombNodePL.h"

using util::geo::Point;
using namespace octi::graph;

// _____________________________________________________________________________
CombNodePL::CombNodePL(octi::graph::Node* parent) : _parent(parent) {}

// _____________________________________________________________________________
const Point* CombNodePL::getGeom() const { return _parent->pl().getGeom(); }

// _____________________________________________________________________________
octi::graph::Node* CombNodePL::getParent() const { return _parent; }

// _____________________________________________________________________________
void CombNodePL::getAttrs(json::object_t& obj) const {
  return _parent->pl().getAttrs(obj);
}

// _____________________________________________________________________________
void CombNodePL::addOrderedEdge(util::graph::Edge<CombNodePL, CombEdgePL>* e,
                                double deg) {
  _edgeOrder.insert(
      std::pair<util::graph::Edge<CombNodePL, CombEdgePL>*, double>(e, deg));
}

// _____________________________________________________________________________
int64_t CombNodePL::distBetween(
    util::graph::Edge<CombNodePL, CombEdgePL>* a,
    util::graph::Edge<CombNodePL, CombEdgePL>* b) const {
  auto aIt = _edgeOrder.begin();
  auto bIt = _edgeOrder.begin();

  for (; aIt != _edgeOrder.end(); aIt++) {
    if (aIt->first == a) break;
  }

  for (; bIt != _edgeOrder.end(); bIt++) {
    if (bIt->first == b) break;
  }

  assert(aIt != _edgeOrder.end());
  assert(bIt != _edgeOrder.end());

  int32_t ap = std::distance(_edgeOrder.begin(), aIt);
  int32_t bp = std::distance(_edgeOrder.begin(), bIt);

  int32_t ret = ((ap > bp ? -1 : 1) * (abs(bp - ap)));

  if (ret < 0) ret = (_edgeOrder.size() + ret);

  return ret;
}

// _____________________________________________________________________________
bool CombNodePL::hasOrderedEdge(
    util::graph::Edge<CombNodePL, CombEdgePL>* e) const {
  auto aIt = _edgeOrder.begin();

  for (; aIt != _edgeOrder.end(); aIt++) {
    if (aIt->first == e) return true;
  }

  return false;
}

// _____________________________________________________________________________
void CombNodePL::setRouteNumber(size_t n) { _routeNumber = n; }

// _____________________________________________________________________________
size_t CombNodePL::getRouteNumber() const { return _routeNumber; }

// _____________________________________________________________________________
std::set<std::pair<util::graph::Edge<CombNodePL, CombEdgePL>*, double>,
               PairCmp>&
CombNodePL::getOrderedEdges() {
  return _edgeOrder;
}

