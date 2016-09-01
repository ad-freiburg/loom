// Copyright 2016, Universityjof Freiburg,
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
    _ordering.push_back(_trips.size() - 1);
  }
  to->addTrip(t, dirNode);
}

// _____________________________________________________________________________
const std::vector<TripOccurance>& EdgeTripGeom::getTripsUnordered()
const {
  return _trips;
}

// _____________________________________________________________________________
TripOccWithPos EdgeTripGeom::getTripsForRoute(const gtfs::Route* r) const {
  for (size_t i = 0; i < _trips.size(); i++) {
    const TripOccurance& to = _trips[i];
    if (to.route == r) {
      size_t pos = std::find(_ordering.begin(), _ordering.end(), i) - _ordering.begin();
      return std::pair<TripOccurance*, size_t>(const_cast<TripOccurance*>(&to), pos);
    }
  }
  return std::pair<TripOccurance*, size_t>(0, 0);
}

// _____________________________________________________________________________
const std::vector<size_t>& EdgeTripGeom::getTripOrdering() const {
  return _ordering;
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
      it = removeTripOccurance(it);
      it--;
    }
  }
}
// _____________________________________________________________________________
std::vector<TripOccurance>::iterator
EdgeTripGeom::removeTripOccurance(std::vector<TripOccurance>::const_iterator pos) {
  std::vector<TripOccurance>::iterator ret = _trips.erase(pos);

  // rebuild ordering
  _ordering.clear();
  for (size_t i = 0; i < _trips.size(); ++i) {
    _ordering.push_back(i);
  }

  return ret;
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

// _____________________________________________________________________________
void EdgeTripGeom::setTripOrdering(std::vector<size_t>& ordering) {
  assert(ordering.size() == _trips.size());
  _ordering = ordering;
}
