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
  if (_is.size() > 0) {
    obj["station_id"] = _is.begin()->id;
    obj["station_label"] = _is.begin()->name;
  }
}

// _____________________________________________________________________________
void NodePL::addStop(StationInfo i) {
  _is.push_back(i);
}

// _____________________________________________________________________________
const std::vector<StationInfo>& NodePL::getStops() const {
  return _is;
}
