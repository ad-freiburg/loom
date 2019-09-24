// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/GridEdgePL.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using namespace octi::gridgraph;

// _____________________________________________________________________________
GridEdgePL::GridEdgePL(double c, bool secondary)
    : _c(c), _isSecondary(secondary), _closed(false), _visited(0) {}

// _____________________________________________________________________________
GridEdgePL::GridEdgePL(double c, bool secondary,
                       bool closed)
    : _c(c), _isSecondary(secondary), _closed(closed), _visited(0) {}

// _____________________________________________________________________________
const util::geo::Line<double>* GridEdgePL::getGeom() const {
  return 0;
}

// _____________________________________________________________________________
void GridEdgePL::addResidentEdge(
    util::graph::Edge<octi::combgraph::CombNodePL, octi::combgraph::CombEdgePL>*
        e) {
  _resEdges.insert(e);
}

// _____________________________________________________________________________
void GridEdgePL::reset() {
  _resEdges.clear();
  _closed = false;
  _visited = 0;
}

// _____________________________________________________________________________
const std::set<util::graph::Edge<octi::combgraph::CombNodePL,
                                 octi::combgraph::CombEdgePL>*>&
GridEdgePL::getResEdges() const {
  return _resEdges;
}

// _____________________________________________________________________________
util::json::Dict GridEdgePL::getAttrs() const {
  util::json::Dict obj;
  std::vector<std::string> routes;
  obj["cost"] = cost() == std::numeric_limits<double>::infinity()
                    ? "inf"
                    : util::toString(cost());
  obj["visited"] = _visited;

  for (auto r : _resEdges) {
    routes.push_back(util::toString(r));
  }
  obj["routes"] = util::implode(routes.begin(), routes.end(), ",");

  return obj;
}

// _____________________________________________________________________________
double GridEdgePL::cost() const {
  return _closed ? std::numeric_limits<double>::infinity() : rawCost();
}

// _____________________________________________________________________________
double GridEdgePL::rawCost() const { return _c; }

// _____________________________________________________________________________
void GridEdgePL::close() { _closed = true; }

// _____________________________________________________________________________
bool GridEdgePL::closed() const { return _closed; }

// _____________________________________________________________________________
void GridEdgePL::open() { _closed = false; }

// _____________________________________________________________________________
void GridEdgePL::setCost(double c) { _c = c; }

// _____________________________________________________________________________
bool GridEdgePL::isSecondary() const { return _isSecondary; }

// _____________________________________________________________________________
void GridEdgePL::setVisited(int i) { _visited = i; }
