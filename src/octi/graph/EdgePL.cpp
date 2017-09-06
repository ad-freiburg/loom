// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/PolyLine.h"
#include "octi/gridgraph/EdgePL.h"

using util::geo::PolyLine;
using namespace octi::gridgraph;

// _____________________________________________________________________________
EdgePL::EdgePL(const PolyLine& pl, double c) : _pl(pl), _c(c) {

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
  obj["routes"] = json::array();
  obj["cost"] = _c;

  for (auto r : _routes) {
    obj["routes"].push_back(r);
  }
}

// _____________________________________________________________________________
double EdgePL::cost() const {
  return _c;
}

// _____________________________________________________________________________
void EdgePL::setCost(double c) {
  _c = c;
}
