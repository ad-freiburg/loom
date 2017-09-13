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

