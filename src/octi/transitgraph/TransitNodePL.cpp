// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/transitgraph/TransitNodePL.h"

using util::geo::Point;
using namespace octi::transitgraph;

// _____________________________________________________________________________
TransitNodePL::TransitNodePL(Point<double> pos) : _pos(pos) {}

// _____________________________________________________________________________
const Point<double>* TransitNodePL::getGeom() const { return &_pos; }

// _____________________________________________________________________________
void TransitNodePL::setGeom(const Point<double>& p) { _pos = p; }

// _____________________________________________________________________________
util::json::Dict TransitNodePL::getAttrs() const {
  util::json::Dict obj;
  if (_is.size() > 0) {
    obj["station_id"] = _is.begin()->id;
    obj["station_label"] = _is.begin()->name;
  }
  return obj;
}

// _____________________________________________________________________________
void TransitNodePL::addStop(StationInfo i) { _is.push_back(i); }

// _____________________________________________________________________________
const std::vector<StationInfo>& TransitNodePL::getStops() const { return _is; }
