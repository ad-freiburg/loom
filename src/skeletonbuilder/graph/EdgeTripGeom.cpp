// Copyright 2016, Universityjof Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "gtfsparser/gtfs/Trip.h"
#include "gtfsparser/gtfs/Route.h"
#include "EdgeTripGeom.h"

using namespace skeletonbuilder;
using namespace graph;
using namespace gtfsparser;

// _____________________________________________________________________________
EdgeTripGeom::EdgeTripGeom(geo::PolyLine geom, const Node* geomDir)
: _geomDir(geomDir) {
  setGeom(geom);
}

// _____________________________________________________________________________
void EdgeTripGeom::addTrip(gtfs::Trip* t, const Node* dirNode,
    geo::PolyLine& pl) {
  std::vector<const geo::PolyLine*> vec;
  vec.push_back(&_geom);
  if (dirNode != _geomDir) {
    pl.reverse();
  }
  vec.push_back(&pl);
  setGeom(geo::PolyLine::average(vec));

  addTrip(t, dirNode);
}

// _____________________________________________________________________________
void EdgeTripGeom::addTrip(gtfs::Trip* t, const Node* dirNode) {
  TripOccurance* to = getTripsForRoute(t->getRoute());
  if (!to) {
    _trips.push_back(TripOccurance(t->getRoute()));
    to = &_trips.back();
  }
  to->addTrip(t, dirNode);
}

// _____________________________________________________________________________
const std::vector<TripOccurance>& EdgeTripGeom::getTripsUnordered()
const {
  return _trips;
}

// _____________________________________________________________________________
std::vector<TripOccurance>* EdgeTripGeom::getTripsUnordered() {
  return &_trips;
}

// _____________________________________________________________________________
TripOccurance* EdgeTripGeom::getTripsForRoute(const gtfs::Route* r) const {
  for (size_t i = 0; i < _trips.size(); i++) {
    TripOccurance* to = const_cast<TripOccurance*>(&_trips[i]);
    if (to->route == r) {
      return to;
    }
  }
  return 0;
}

// _____________________________________________________________________________
const geo::PolyLine& EdgeTripGeom::getGeom() const {
  return _geom;
}

// _____________________________________________________________________________
void EdgeTripGeom::setGeom(const geo::PolyLine& p) {
  _geom = p;
}

// _____________________________________________________________________________
bool EdgeTripGeom::containsRoute(gtfs::Route* r) const {
  if (getTripsForRoute(r)) return true;

  return false;
}

// _____________________________________________________________________________
size_t EdgeTripGeom::getTripCardinality() const {
  size_t ret = 0;

  for (auto& t : _trips) {
    ret += t.trips.size();
  }

  return ret;
}

// _____________________________________________________________________________
size_t EdgeTripGeom::getCardinality() const {
  return _trips.size();
}

// _____________________________________________________________________________
void EdgeTripGeom::removeOrphans() {
  double avgTrips = 0;

  for (auto& to : _trips) {
    avgTrips += to.trips.size();
  }
  avgTrips /= _trips.size();

  for (auto it = _trips.begin(); it != _trips.end(); it++) {
    if (it->trips.size() < avgTrips * 0.15) {
      it = removeTripOccurance(it);
      it--;
    }
  }
}
// _____________________________________________________________________________
std::vector<TripOccurance>::iterator
EdgeTripGeom::removeTripOccurance(std::vector<TripOccurance>::const_iterator pos) {
  return _trips.erase(pos);
}

// _____________________________________________________________________________
const Node* EdgeTripGeom::getGeomDir() const {
  return _geomDir;
}

// _____________________________________________________________________________
void EdgeTripGeom::setGeomDir(const Node* n) {
  _geomDir = n;
}

// _____________________________________________________________________________
bool EdgeTripGeom::routeEquivalent(const EdgeTripGeom& g) const {
  if (_trips.size() != g._trips.size()) return false;

  for (auto to : _trips) {
    if (!g.containsRoute(to.route)) {
      return false;
    }
  }

  return true;
}