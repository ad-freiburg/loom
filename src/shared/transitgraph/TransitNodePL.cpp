// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/transitgraph/TransitNodePL.h"

using util::geo::Point;
using namespace shared::transitgraph;

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

  auto arr = util::json::Array();

  for (const auto& ro : _connEx) {
    for (const auto& exFr : ro.second) {
      for (const auto* exTo : exFr.second) {
        obj["route"] = util::toString(ro.first->getId());
        obj["edge1_node"] = util::toString(exFr.first);
        obj["edge2_node"] = util::toString(exTo);
        arr.push_back(obj);
      }
    }
  }

  if (arr.size()) obj["excluded_line_conns"] = arr;
  return obj;
}

// _____________________________________________________________________________
void TransitNodePL::addStop(StationInfo i) { _is.push_back(i); }

// _____________________________________________________________________________
const std::vector<StationInfo>& TransitNodePL::getStops() const { return _is; }

// _____________________________________________________________________________
void TransitNodePL::addConnExc(const Route* r, const TransitEdge* edgeA,
                               const TransitEdge* edgeB) {
  _connEx[r][edgeA].insert(edgeB);
  // index the other direction also, will lead to faster lookups later on
  _connEx[r][edgeB].insert(edgeA);
}

// _____________________________________________________________________________
bool TransitNodePL::connOccurs(const Route* r, const TransitEdge* edgeA,
                               const TransitEdge* edgeB) const {
  const auto& i = _connEx.find(r);
  if (_connEx.find(r) == _connEx.end()) return true;

  const auto& ii = i->second.find(edgeA);
  if (ii == i->second.end()) return true;

  return ii->second.find(edgeB) == ii->second.end();
}
