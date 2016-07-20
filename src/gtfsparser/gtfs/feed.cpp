// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <unordered_map>
#include <string>
#include "feed.h"
#include "agency.h"
#include "route.h"
#include "trip.h"
#include "shape.h"

using namespace gtfsparser;
using namespace gtfs;

// ____________________________________________________________________________
bool Feed::addAgency(Agency* a) {
  return (_agencies.insert(std::pair<std::string, Agency*>(a->getId(), a))).second;
}

// ____________________________________________________________________________
const Agency* Feed::getAgencyById(const std::string& id) const {
  if (_agencies.find(id) != _agencies.end()) {
    return _agencies.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
Agency* Feed::getAgencyById(const std::string& id) {
  if (_agencies.find(id) != _agencies.end()) {
    return _agencies.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
bool Feed::addStop(Stop* a) {
  return (_stops.insert(std::pair<std::string, Stop*>(a->getId(), a))).second;
}

// ____________________________________________________________________________
const Stop* Feed::getStopById(const std::string& id) const {
  if (_stops.find(id) != _stops.end()) {
    return _stops.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
Stop* Feed::getStopById(const std::string& id) {
  if (_stops.find(id) != _stops.end()) {
    return _stops.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
bool Feed::addRoute(Route* r) {
  return (_routes.insert(std::pair<std::string, Route*>(r->getId(), r))).second;
}

// ____________________________________________________________________________
const Route* Feed::getRouteById(const std::string& id) const {
  if (_routes.find(id) != _routes.end()) {
    return _routes.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
Route* Feed::getRouteById(const std::string& id) {
  if (_routes.find(id) != _routes.end()) {
    return _routes.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
bool Feed::addTrip(Trip* t) {
  return (_trips.insert(std::pair<std::string, Trip*>(t->getId(), t))).second;
}

// ____________________________________________________________________________
const Trip* Feed::getTripById(const std::string& id) const {
  if (_trips.find(id) != _trips.end()) {
    return _trips.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
Trip* Feed::getTripById(const std::string& id) {
  if (_trips.find(id) != _trips.end()) {
    return _trips.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
bool Feed::addShape(Shape* s) {
  return (_shapes.insert(std::pair<std::string, Shape*>(s->getId(), s))).second;
}

// ____________________________________________________________________________
const Shape* Feed::getShapeById(const std::string& id) const {
  if (_shapes.find(id) != _shapes.end()) {
    return _shapes.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
Shape* Feed::getShapeById(const std::string& id) {
  if (_shapes.find(id) != _shapes.end()) {
    return _shapes.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
bool Feed::addService(Service* s) {
  return (_services.insert(std::pair<std::string, Service*>(s->getId(), s))).second;
}

// ____________________________________________________________________________
const Service* Feed::getServiceById(const std::string& id) const {
  if (_services.find(id) != _services.end()) {
    return _services.find(id)->second;
  }
  return 0;
}

// ____________________________________________________________________________
Service* Feed::getServiceById(const std::string& id) {
  if (_services.find(id) != _services.end()) {
    return _services.find(id)->second;
  }
  return 0;
}

