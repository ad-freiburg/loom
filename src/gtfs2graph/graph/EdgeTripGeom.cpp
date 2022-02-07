// Copyright 2016, Universityjof Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "ad/cppgtfs/gtfs/Route.h"
#include "ad/cppgtfs/gtfs/Trip.h"
#include "gtfs2graph/graph/EdgeTripGeom.h"

using namespace gtfs2graph;
using namespace graph;
using namespace ad::cppgtfs;
using namespace util::geo;

// _____________________________________________________________________________
EdgeTripGeom::EdgeTripGeom(PolyLine<double> geom, const Node* geomDir)
    : _geomDir(geomDir) {
  setGeom(geom);
}

// _____________________________________________________________________________
void EdgeTripGeom::addTrip(gtfs::Trip* t, const Node* dirNode,
                           PolyLine<double>& pl) {
  if (dirNode != _geomDir) pl.reverse();

  // TODO: only do this if they are equal within some threshold!
  setGeom(PolyLine<double>::average({&_geom, &pl}));

  addTrip(t, dirNode);
}

// _____________________________________________________________________________
void EdgeTripGeom::addTrip(gtfs::Trip* t, const Node* dirNode) {
  RouteOccurance* to = getRouteOcc(t->getRoute());
  if (!to) {
    _routeOccs.push_back(RouteOccurance(t->getRoute()));
    to = &_routeOccs.back();
  }
  to->addTrip(t, dirNode);
}

// _____________________________________________________________________________
const std::vector<RouteOccurance>& EdgeTripGeom::getTripsUnordered() const {
  return _routeOccs;
}

// _____________________________________________________________________________
std::vector<RouteOccurance>* EdgeTripGeom::getTripsUnordered() {
  return &_routeOccs;
}

// _____________________________________________________________________________
RouteOccurance* EdgeTripGeom::getRouteOcc(const gtfs::Route* r) const {
  for (size_t i = 0; i < _routeOccs.size(); i++) {
    RouteOccurance* to = const_cast<RouteOccurance*>(&_routeOccs[i]);
    if (to->route == r) return to;
  }
  return 0;
}

// _____________________________________________________________________________
const PolyLine<double>& EdgeTripGeom::getGeom() const { return _geom; }

// _____________________________________________________________________________
void EdgeTripGeom::setGeom(const PolyLine<double>& p) { _geom = p; }

// _____________________________________________________________________________
bool EdgeTripGeom::containsRoute(gtfs::Route* r) const {
  return getRouteOcc(r);
}

// _____________________________________________________________________________
size_t EdgeTripGeom::getTripCardinality() const {
  size_t ret = 0;

  for (auto& t : _routeOccs) ret += t.trips.size();

  return ret;
}

// _____________________________________________________________________________
size_t EdgeTripGeom::getCardinality() const { return _routeOccs.size(); }

// _____________________________________________________________________________
std::vector<RouteOccurance>::iterator EdgeTripGeom::removeRouteOcc(
    std::vector<RouteOccurance>::const_iterator pos) {
  return _routeOccs.erase(pos);
}

// _____________________________________________________________________________
const Node* EdgeTripGeom::getGeomDir() const { return _geomDir; }

// _____________________________________________________________________________
void EdgeTripGeom::setGeomDir(const Node* n) { _geomDir = n; }
