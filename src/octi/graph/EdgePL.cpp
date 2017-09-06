// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/PolyLine.h"
#include "octi/gridgraph/EdgePL.h"

using util::geo::PolyLine;
using namespace octi::gridgraph;

// _____________________________________________________________________________
EdgePL::EdgePL() {

}

// _____________________________________________________________________________
const util::geo::Line* EdgePL::getGeom() const {
  return &_pl.getLine();
}

// _____________________________________________________________________________
void EdgePL::addRoute(std::string r) {
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