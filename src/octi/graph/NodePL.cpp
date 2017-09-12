// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/NodePL.h"

using util::geo::Point;
using namespace octi::graph;

// _____________________________________________________________________________
NodePL::NodePL(Point pos) : _pos(pos) {
}

// _____________________________________________________________________________
const Point* NodePL::getGeom() const {
  return &_pos;
}

// _____________________________________________________________________________
void NodePL::setGeom(const Point& p) {
   _pos = p;
}

// _____________________________________________________________________________
void NodePL::getAttrs(json::object_t& obj) const {

}

