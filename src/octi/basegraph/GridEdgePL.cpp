// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/basegraph/GridEdgePL.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using namespace octi::basegraph;

// _____________________________________________________________________________
GridEdgePL::GridEdgePL(double c, bool secondary)
    : _c(c),
      _isSecondary(secondary),
      _closed(false),
      _blocked(false),
      _resEdgs(0), _rndrOrder(0) {}

// _____________________________________________________________________________
GridEdgePL::GridEdgePL(double c, bool secondary, bool closed)
    : _c(c),
      _isSecondary(secondary),
      _closed(closed),
      _blocked(false),
      _resEdgs(0), _rndrOrder(0) {}

// _____________________________________________________________________________
const util::geo::Line<double>* GridEdgePL::getGeom() const { return 0; }

// _____________________________________________________________________________
void GridEdgePL::reset() {
  _closed = false;
  clearResEdges();
}

// _____________________________________________________________________________
util::json::Dict GridEdgePL::getAttrs() const {
  util::json::Dict obj;
  obj["cost"] = cost() == std::numeric_limits<double>::infinity()
                    ? "inf"
                    : util::toString(cost());
  obj["res_edges"] = util::toString((int)_resEdgs);
  obj["rndr_order"] = util::toString((int)_rndrOrder);
  return obj;
}
// _____________________________________________________________________________
double GridEdgePL::cost() const {
  return (_closed || _blocked) ? std::numeric_limits<double>::infinity()
                               : rawCost();
}

// _____________________________________________________________________________
double GridEdgePL::rawCost() const { return _c; }

// _____________________________________________________________________________
void GridEdgePL::addResEdge() { _resEdgs++; }

// _____________________________________________________________________________
void GridEdgePL::close() { _closed = true; }

// _____________________________________________________________________________
bool GridEdgePL::closed() const { return _closed; }

// _____________________________________________________________________________
void GridEdgePL::open() { _closed = false; }

// _____________________________________________________________________________
void GridEdgePL::block() { _blocked = true; }

// _____________________________________________________________________________
void GridEdgePL::unblock() { _blocked = false; }

// _____________________________________________________________________________
void GridEdgePL::setCost(double c) { _c = c; }

// _____________________________________________________________________________
bool GridEdgePL::isSecondary() const { return _isSecondary; }

// _____________________________________________________________________________
void GridEdgePL::clearResEdges() { _resEdgs = 0; }

// _____________________________________________________________________________
void GridEdgePL::setId(size_t id) { _id = id; }

// _____________________________________________________________________________
size_t GridEdgePL::getId() const { return _id; }

// _____________________________________________________________________________
void GridEdgePL::setRndrOrder(size_t order) {
  _rndrOrder = order;
}
