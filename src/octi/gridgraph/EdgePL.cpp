// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/EdgePL.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using namespace octi::gridgraph;

// _____________________________________________________________________________
EdgePL::EdgePL(const PolyLine<double>& pl, double c, bool secondary)
    : _pl(pl), _c(c), _isSecondary(secondary), _closed(false), _visited(0) {}

// _____________________________________________________________________________
EdgePL::EdgePL(const PolyLine<double>& pl, double c, bool secondary, bool closed)
    : _pl(pl), _c(c), _isSecondary(secondary), _closed(closed), _visited(0) {}

// _____________________________________________________________________________
const util::geo::Line<double>* EdgePL::getGeom() const { return &_pl.getLine(); }

// _____________________________________________________________________________
void EdgePL::addResidentEdge(
    util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>* e) {
  _resEdges.insert(e);
}
//
// _____________________________________________________________________________
const std::set<
    util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>*>&
EdgePL::getResEdges() const {
  return _resEdges;
}

// _____________________________________________________________________________
util::json::Dict EdgePL::getAttrs() const {
  util::json::Dict obj;
  util::json::Array routes;
  obj["cost"] = cost() == std::numeric_limits<double>::infinity()
                    ? "inf"
                    : util::toString(cost());
  obj["visited"] = _visited;

  for (auto r : _resEdges) {
    routes.push_back(util::toString(r));
  }
  obj["routes"] = routes;

  return obj;
}

// _____________________________________________________________________________
double EdgePL::cost() const {
  return _closed ? std::numeric_limits<double>::infinity() : rawCost();
}

// _____________________________________________________________________________
double EdgePL::rawCost() const { return _c; }

// _____________________________________________________________________________
void EdgePL::close() { _closed = true; }

// _____________________________________________________________________________
bool EdgePL::closed() const { return _closed; }

// _____________________________________________________________________________
void EdgePL::open() { _closed = false; }

// _____________________________________________________________________________
void EdgePL::setCost(double c) { _c = c; }

// _____________________________________________________________________________
bool EdgePL::isSecondary() const { return _isSecondary; }

// _____________________________________________________________________________
void EdgePL::setVisited(int i) { _visited = i; }
