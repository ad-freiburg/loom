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
void EdgeTripGeom::addTrip(gtfs::Trip* t) {
  _trips[t->getRoute()].push_back(t);
}

// _____________________________________________________________________________
const std::map<gtfs::Route*, std::vector<gtfs::Trip*> >& EdgeTripGeom::getTrips() 
const {
  return _trips;
}

// _____________________________________________________________________________
const geo::PolyLine& EdgeTripGeom::getGeom() const {
  return _geom;
}
