// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/CombNodePL.h"
#include "octi/graph/EdgeOrdering.h"
using util::geo::Point;
using namespace octi::graph;

// _____________________________________________________________________________
CombNodePL::CombNodePL(octi::graph::Node* parent) : _parent(parent) {}

// _____________________________________________________________________________
const Point<double>* CombNodePL::getGeom() const {
  return _parent->pl().getGeom();
}

// _____________________________________________________________________________
octi::graph::Node* CombNodePL::getParent() const { return _parent; }

// _____________________________________________________________________________
util::json::Dict CombNodePL::getAttrs() const {
  return _parent->pl().getAttrs();
}

// _____________________________________________________________________________
void CombNodePL::setRouteNumber(size_t n) { _routeNumber = n; }

// _____________________________________________________________________________
size_t CombNodePL::getRouteNumber() const { return _routeNumber; }

// _____________________________________________________________________________
const EdgeOrdering& CombNodePL::getEdgeOrdering() { return _ordering; }

// _____________________________________________________________________________
void CombNodePL::setEdgeOrdering(const EdgeOrdering& ordering) {
  _ordering = ordering;
}

// _____________________________________________________________________________
std::string CombNodePL::toString() const {
  std::stringstream ret;
  ret << "<" << this << " (" << getParent() << ")";
  if (getParent()->pl().getStops().size()) {
    ret << " " << getParent()->pl().getStops().front().name;
  }

  ret << ">";

  return ret.str();
}
