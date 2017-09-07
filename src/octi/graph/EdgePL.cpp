// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/PolyLine.h"
#include "octi/graph/EdgePL.h"

using util::geo::PolyLine;
using namespace octi::graph;

// _____________________________________________________________________________
EdgePL::EdgePL(const PolyLine& p) : _p(p) {

}

// _____________________________________________________________________________
const util::geo::Line* EdgePL::getGeom() const {
  return &_p.getLine();
}

// _____________________________________________________________________________
void EdgePL::addRoute(const std::string& r) {
  _routes.insert(r);
}
//
// _____________________________________________________________________________
const std::set<std::string>& EdgePL::getRoutes() const {
  return _routes;
}

// _____________________________________________________________________________
void EdgePL::getAttrs(json::object_t& obj) const {

}
