// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "gtfsparser/gtfs/Trip.h"
#include "gtfsparser/gtfs/Route.h"
#include "EdgeTripGeom.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;

// _____________________________________________________________________________
EdgeTripGeom::EdgeTripGeom(geo::PolyLine geom, const Node* geomDir)
: _w(1), _s(1), _geomDir(geomDir) {
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
  _geom = geo::PolyLine::average(vec);

  addTrip(t, dirNode);
}

// _____________________________________________________________________________
void EdgeTripGeom::addTrip(gtfs::Trip* t, const Node* dirNode) {
  TripOccurance* to = getTripsForRoute(t->getRoute()).first;
  if (!to) {
    _trips.push_back(TripOccurance(t->getRoute()));
    to = &_trips.back();
  }
  to->addTrip(t, dirNode);
}

// _____________________________________________________________________________
const std::vector<TripOccurance>& EdgeTripGeom::getTrips()
const {
  return _trips;
}

// _____________________________________________________________________________
TripOccWithPos EdgeTripGeom::getTripsForRoute(const gtfs::Route* r) const {
  for (size_t i = 0; i < _trips.size(); i++) {
    const TripOccurance& to = _trips[i];
    if (to.route == r)
      return std::pair<TripOccurance*, size_t>(const_cast<TripOccurance*>(&to), i);
  }
  return std::pair<TripOccurance*, size_t>(0, 0);
}

// _____________________________________________________________________________
std::vector<TripOccurance>* EdgeTripGeom::getTrips() {
  return &_trips;
}

// _____________________________________________________________________________
const geo::PolyLine& EdgeTripGeom::getGeom() const {
  return _geom;
}

// _____________________________________________________________________________
void EdgeTripGeom::setGeom(const geo::PolyLine& p) {
  _geom = p;
  if (util::geo::dist(_geom.getLine().back(), _geomDir->getPos()) <= util::geo::dist(_geom.getLine().front(), _geomDir->getPos())) {
 //   _geom.reverse();
  }
}

// _____________________________________________________________________________
bool EdgeTripGeom::containsRoute(gtfs::Route* r) const {
  return getTripsForRoute(r).first;
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
void EdgeTripGeom::removeOrphans() {
  double avgTrips = 0;

  for (auto& to : _trips) {
    avgTrips += to.trips.size();
  }
  avgTrips /= _trips.size();

  for (auto it = _trips.begin(); it != _trips.end(); it++) {
    if (it->trips.size() < avgTrips * 0.1) {
      it = _trips.erase(it);
      it--;
    }
  }
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
double EdgeTripGeom::getWidth() const {
  return _w;
}

// _____________________________________________________________________________
double EdgeTripGeom::getSpacing() const {
  return _s;
}

// _____________________________________________________________________________
void EdgeTripGeom::setWidth(double w) {
  _w = w;
}

// _____________________________________________________________________________
void EdgeTripGeom::setSpacing(double s) {
  _s = s;
}

// _____________________________________________________________________________
double EdgeTripGeom::getTotalWidth() const {
  return _w * _trips.size() + _s * (_trips.size() - 1);
}
