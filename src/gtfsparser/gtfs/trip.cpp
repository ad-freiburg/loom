// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <iostream>
#include <stdint.h>
#include <string>
#include "route.h"
#include "service.h"
#include "shape.h"
#include "stoptime.h"
#include "trip.h"

using std::exception;
using std::string;


using namespace gtfsparser;
using namespace gtfs;

// _____________________________________________________________________________
Trip::Trip(std::string id, Route* r, Service* s, std::string hs,
  std::string short_name, DIRECTION dir, std::string blockid, Shape* shp,
  WC_BIKE_ACCESSIBLE wc, WC_BIKE_ACCESSIBLE ba)
: _id(id), _route(r), _service(s), _headsign(hs), _short_name(short_name),
  _dir(dir), _block_id(blockid), _shape(shp), _wc(wc), _ba(ba) {}

// _____________________________________________________________________________
const std::string& Trip::getId() const {
  return _id;
}

// _____________________________________________________________________________
Route* Trip::getRoute() const {
  return _route;
}

// _____________________________________________________________________________
Service* Trip::getService() const {
  return _service;
}

// _____________________________________________________________________________
const std::string& Trip::getHeadsign() const {
  return _headsign;
}

// _____________________________________________________________________________
const std::string& Trip::getShortname() const {
  return _short_name;
}

// _____________________________________________________________________________
Trip::DIRECTION Trip::getDirection() const {
  return _dir;
}

// _____________________________________________________________________________
const std::string& Trip::getBlockId() const {
  return _block_id;
}

// _____________________________________________________________________________
Shape* Trip::getShape() const {
  return _shape;
}

// _____________________________________________________________________________
Trip::WC_BIKE_ACCESSIBLE Trip::getWheelchairAccessibility() const {
  return _wc;
}

// _____________________________________________________________________________
Trip::WC_BIKE_ACCESSIBLE Trip::getBikesAllowed() const {
  return _ba;
}

// _____________________________________________________________________________
const StopTimes& Trip::getStopTimes() const {
  return _stoptimes;
}

// _____________________________________________________________________________
bool Trip::addStopTime(const StopTime& t) {
  return _stoptimes.insert(t).second ;
}
