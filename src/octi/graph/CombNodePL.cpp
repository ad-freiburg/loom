// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/CombNodePL.h"

using util::geo::Point;
using namespace octi::graph;

// _____________________________________________________________________________
CombNodePL::CombNodePL(octi::graph::Node* parent) : _parent(parent) {
}

// _____________________________________________________________________________
const Point* CombNodePL::getGeom() const {
  return _parent->pl().getGeom();
}

// _____________________________________________________________________________
octi::graph::Node* CombNodePL::getParent() const {
  return _parent;
}

// _____________________________________________________________________________
void CombNodePL::getAttrs(json::object_t& obj) const {
  return _parent->pl().getAttrs(obj);
}

// _____________________________________________________________________________
void CombNodePL::addOrderedEdge(util::graph::Edge<CombNodePL, CombEdgePL>* e, double deg) {
  _edgeOrder.insert(std::pair<util::graph::Edge<CombNodePL, CombEdgePL>*, double>(e, deg));
}

// _____________________________________________________________________________
int32_t CombNodePL::distBetween(util::graph::Edge<CombNodePL, CombEdgePL>* a, util::graph::Edge<CombNodePL, CombEdgePL>* b) const {
  return 0;
}

