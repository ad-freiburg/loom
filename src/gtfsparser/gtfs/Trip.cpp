// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>
#include <iostream>
#include <string>
#include "Route.h"
#include "Service.h"
#include "Shape.h"
#include "StopTime.h"
#include "Trip.h"

using std::exception;
using std::string;

using gtfsparser::gtfs::Trip;
using gtfsparser::gtfs::Route;
using gtfsparser::gtfs::Service;
using gtfsparser::gtfs::Shape;
using gtfsparser::gtfs::StopTime;
using gtfsparser::gtfs::StopTimes;

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
  return _stoptimes.insert(t).second;
}
