// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/NodePL.h"

using util::geo::Point;
using namespace octi::graph;

// _____________________________________________________________________________
NodePL::NodePL(Point<double> pos) : _pos(pos) {
}

// _____________________________________________________________________________
const Point<double>* NodePL::getGeom() const {
  return &_pos;
}

// _____________________________________________________________________________
void NodePL::setGeom(const Point<double>& p) {
   _pos = p;
}

// _____________________________________________________________________________
util::json::Dict NodePL::getAttrs() const {
  util::json::Dict obj;
  if (_is.size() > 0) {
    obj["station_id"] = _is.begin()->id;
    obj["station_label"] = _is.begin()->name;
  }
  return obj;
}

// _____________________________________________________________________________
void NodePL::addStop(StationInfo i) {
  _is.push_back(i);
}

// _____________________________________________________________________________
const std::vector<StationInfo>& NodePL::getStops() const {
  return _is;
}
