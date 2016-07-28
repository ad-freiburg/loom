// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <vector>
#include "gtfsparser/gtfs/trip.h"
#include "gtfsparser/gtfs/route.h"
#include "edgetripgeom.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;

// _____________________________________________________________________________
EdgeTripGeom::EdgeTripGeom(geo::PolyLine geom) : _geom(geom) {

}

// _____________________________________________________________________________
void EdgeTripGeom::addTrip(gtfs::Trip* t, Node* dirNode) {
  _trips[t->getRoute()].addTrip(t, dirNode);
}

// _____________________________________________________________________________
const std::map<gtfs::Route*, TripOccurance>& EdgeTripGeom::getTrips()
const {
  return _trips;
}

// _____________________________________________________________________________
std::map<gtfs::Route*, TripOccurance >* EdgeTripGeom::getTrips() {
  return &_trips;
}

// _____________________________________________________________________________
const geo::PolyLine& EdgeTripGeom::getGeom() const {
  return _geom;
}

// _____________________________________________________________________________
bool EdgeTripGeom::containsRoute(gtfs::Route* r) const {
  return _trips.find(r) != _trips.end();
}

// _____________________________________________________________________________
size_t EdgeTripGeom::getTripCardinality() const {
  size_t ret = 0;

  for (auto& t : _trips) {
    ret += t.second.trips.size();
  }

  return ret;
}

// _____________________________________________________________________________
void EdgeTripGeom::removeOrphans() {
  double avgTrips = 0;

  for (auto& to : _trips) {
    avgTrips += to.second.trips.size();
  }
  avgTrips /= _trips.size();

  for (auto it = _trips.begin(); it != _trips.end(); it++) {
    if (it->second.trips.size() < avgTrips * 0.1) {
      it = _trips.erase(it);
    }
  }
}
