// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/combgraph/CombNodePL.h"
#include "octi/combgraph/EdgeOrdering.h"
using util::geo::Point;
using octi::combgraph::CombNodePL;

// _____________________________________________________________________________
CombNodePL::CombNodePL(shared::linegraph::LineNode* parent)
    : _parent(parent) {}

// _____________________________________________________________________________
const Point<double>* CombNodePL::getGeom() const {
  return _parent->pl().getGeom();
}

// _____________________________________________________________________________
shared::linegraph::LineNode* CombNodePL::getParent() const {
  return _parent;
}

// _____________________________________________________________________________
util::json::Dict CombNodePL::getAttrs() const {
  return _parent->pl().getAttrs();
}

// _____________________________________________________________________________
void CombNodePL::setRouteNumber(size_t n) { _routeNumber = n; }

// _____________________________________________________________________________
size_t CombNodePL::getLDeg() const { return _routeNumber; }

// _____________________________________________________________________________
const octi::combgraph::EdgeOrdering& CombNodePL::getEdgeOrdering() {
  return _ordering;
}

// _____________________________________________________________________________
void CombNodePL::setEdgeOrdering(
    const octi::combgraph::EdgeOrdering& ordering) {
  _ordering = ordering;
}

// _____________________________________________________________________________
std::string CombNodePL::toString() const {
  std::stringstream ret;
  ret << "<" << this << " (" << getParent() << ")";
  if (getParent()->pl().stops().size()) {
    ret << " " << getParent()->pl().stops().front().name;
  }

  ret << ">";

  return ret.str();
}
