// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/PolyLine.h"
#include "octi/gridgraph/EdgePL.h"
#include "util/String.h"

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
void EdgePL::addResidentEdge(util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>* e) {
  _resEdges.insert(e);
}
//
// _____________________________________________________________________________
const std::set<util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>*>& EdgePL::getResEdges() const {
  return _resEdges;
}

// _____________________________________________________________________________
void EdgePL::getAttrs(json::object_t& obj) const {
  obj["routes"] = json::array();
  obj["cost"] = _c;

  for (auto r : _resEdges) {
    obj["routes"].push_back(util::toString(r));
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
