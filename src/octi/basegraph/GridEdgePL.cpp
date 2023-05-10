// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/basegraph/BaseGraph.h"
#include "octi/basegraph/GridEdgePL.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using namespace octi::basegraph;

// _____________________________________________________________________________
GridEdgePL::GridEdgePL(double c, bool secondary, bool sink)
    : _c(c),
      _isSecondary(secondary),
      _isSink(sink),
      _closed(false),
      _softClosed(false),
      _blocked(false),
      _resEdgs(0) {}

// _____________________________________________________________________________
const util::geo::Line<double>* GridEdgePL::getGeom() const { return 0; }

// _____________________________________________________________________________
size_t GridEdgePL::resEdgs() const { return _resEdgs; }

// _____________________________________________________________________________
void GridEdgePL::reset() {
  _closed = false;
  _resEdgs = 0;
}

// _____________________________________________________________________________
util::json::Dict GridEdgePL::getAttrs() const {
  util::json::Dict obj;
  obj["cost"] = cost() == std::numeric_limits<double>::infinity()
                    ? "inf"
                    : util::toString(cost());
  obj["res_edges"] = util::toString((int)_resEdgs);
  obj["secondary"] = util::toString((int)_isSecondary);
  obj["sink"] = util::toString((int)_isSink);
  obj["closed"] = util::toString(_closed);
  obj["blocked"] = util::toString(_blocked);
  obj["softclosed"] = util::toString(_softClosed);
  return obj;
}
// _____________________________________________________________________________
double GridEdgePL::cost() const {
  // testing relaxed constraints for diagonal intersections
  if (_softClosed || _blocked) return SOFT_INF + rawCost();
  if (_closed) return INF;

  return rawCost();
}

// _____________________________________________________________________________
double GridEdgePL::rawCost() const { return _c; }

// _____________________________________________________________________________
void GridEdgePL::addResEdge() { _resEdgs++; }

// _____________________________________________________________________________
void GridEdgePL::close() {
  _closed = true;
  _softClosed = false;
}

// _____________________________________________________________________________
void GridEdgePL::softClose() {
  if (!_closed) _softClosed = true;
  _closed = true;
}

// _____________________________________________________________________________
bool GridEdgePL::closed() const { return _closed; }

// _____________________________________________________________________________
void GridEdgePL::open() {
  _closed = false;
  _softClosed = false;
}

// _____________________________________________________________________________
void GridEdgePL::block() { _blocked = true; }

// _____________________________________________________________________________
void GridEdgePL::unblock() { _blocked = false; }

// _____________________________________________________________________________
void GridEdgePL::setCost(double c) { _c = c; }

// _____________________________________________________________________________
bool GridEdgePL::isSecondary() const { return _isSecondary; }

// _____________________________________________________________________________
void GridEdgePL::delResEdg() {
  if (_resEdgs > 0) _resEdgs--;
}

// _____________________________________________________________________________
void GridEdgePL::setId(size_t id) { _id = id; }

// _____________________________________________________________________________
size_t GridEdgePL::getId() const { return _id; }
