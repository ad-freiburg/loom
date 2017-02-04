// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "Edge.h"
#include "Node.h"
#include "gtfsparser/gtfs/Trip.h"
#include "EdgeTripGeom.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;

// _____________________________________________________________________________
Edge::Edge(Node* from, Node* to, geo::PolyLine pl, double w,
    double s) : _from(from), _to(to), _width(w), _spacing(s) {
  setGeom(pl);
}

// _____________________________________________________________________________
Node* Edge::getFrom() const {
  return _from;
}

// _____________________________________________________________________________
Node* Edge::getTo() const {
  return _to;
}

// _____________________________________________________________________________
bool Edge::addTrip(gtfs::Trip* t, Node* toNode) {
  assert(toNode == _from || toNode == _to);
  for (auto& e : _tripsContained) {
    if (e.containsRoute(t->getRoute())) {
      for (auto& tr : e.getTripsForRoute(t->getRoute())->trips) {
        // shorcut: if a trip is contained here with the same shape id,
        // don't require recalc of polyline etc
        if (tr->getShape() == t->getShape()) {
          e.addTrip(t, toNode);
          return true;
       }
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
bool Edge::addTrip(gtfs::Trip* t, geo::PolyLine pl, Node* toNode, double w,
    double s) {
  assert(toNode == _from || toNode == _to);

  for (auto& e : _tripsContained) {
    e.addTrip(t, toNode, pl);
    break;
  }

  return true;
}

// _____________________________________________________________________________
const std::vector<TripOccurance>& Edge::getTripsUnordered()
const {
  return _trips;
}

// _____________________________________________________________________________
std::vector<TripOccurance>* Edge::getTripsUnordered() {
  return &_trips;
}

// _____________________________________________________________________________
TripOccurance* Edge::getTripsForRoute(const gtfs::Route* r) const {
  for (size_t i = 0; i < _trips.size(); i++) {
    TripOccurance* to = const_cast<TripOccurance*>(&_trips[i]);
    if (to->route == r) {
      return to;
    }
  }
  return 0;
}

// _____________________________________________________________________________
TripOccWithPos Edge::getTripsForRouteUnder(const gtfs::Route* r,
    const std::vector<size_t> ordering) const {
  for (size_t i = 0; i < _trips.size(); i++) {
    const TripOccurance& to = _trips[i];
    if (to.route == r) {
      size_t pos = std::find(ordering.begin(), ordering.end(), i) - ordering.begin();
      return std::pair<TripOccurance*, size_t>(const_cast<TripOccurance*>(&to), pos);
    }
  }
  return std::pair<TripOccurance*, size_t>(0, 0);
}

// _____________________________________________________________________________
bool Edge::containsRoute(gtfs::Route* r) const {
  if (getTripsForRoute(r)) return true;

  return false;
}

// _____________________________________________________________________________
size_t Edge::getTripCardinality() const {
  size_t ret = 0;

  for (auto& t : _trips) {
    ret += t.trips.size();
  }

  return ret;
}

// _____________________________________________________________________________
size_t Edge::getCardinality() const {
  return _trips.size();
}

// _____________________________________________________________________________
std::vector<TripOccurance>::iterator
Edge::removeTripOccurance(std::vector<TripOccurance>::const_iterator pos) {
  return _trips.erase(pos);
}

// _____________________________________________________________________________
const Node* Edge::getGeomDir() const {
  return _geomDir;
}

// _____________________________________________________________________________
void Edge::setGeomDir(const Node* n) {
  _geomDir = n;
}

// _____________________________________________________________________________
double Edge::getWidth() const {
  return _width;
}

// _____________________________________________________________________________
double Edge::getSpacing() const {
  return _spacing;
}

// _____________________________________________________________________________
double Edge::getTotalWidth() const {
  return getWidth() * _trips.size() + getSpacing() * (_trips.size() - 1);
}

// _____________________________________________________________________________
std::vector<gtfs::Route*> Edge::getSharedRoutes(const Edge& e)
const {
  std::vector<gtfs::Route*> ret;
  for (auto& to : _trips) {
    if (e.containsRoute(to.route)) ret.push_back(to.route);
  }

  return ret;
}

// _____________________________________________________________________________
const geo::PolyLine& Edge::getGeom() const {
  return _geom;
}

// _____________________________________________________________________________
void Edge::setGeom(const geo::PolyLine& p) {
  _geom = p;
}