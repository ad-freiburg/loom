// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include "feed.h"
#include "agency.h"
#include "route.h"

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


