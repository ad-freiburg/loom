// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_FEED_H_
#define GTFSPARSER_GTFS_FEED_H_

#include <unordered_map>
#include <string>
#include "agency.h"
#include "stop.h"
#include "route.h"
#include "trip.h"
#include "shape.h"
#include "service.h"

namespace gtfsparser {
namespace gtfs {

class Feed {
 public:
  Feed() {};

  bool addAgency(Agency* a);
  const Agency* getAgencyById(const std::string& id) const;
  Agency* getAgencyById(const std::string& id);

  bool addStop(Stop* a);
  const Stop* getStopById(const std::string& id) const;
  Stop* getStopById(const std::string& id);

  bool addRoute(Route* a);
  const Route* getRouteById(const std::string& id) const;
  Route* getRouteById(const std::string& id);

  bool addTrip(Trip* a);
  const Trip* getTripById(const std::string& id) const;
  Trip* getTripById(const std::string& id);

  bool addShape(Shape* a);
  const Shape* getShapeById(const std::string& id) const;
  Shape* getShapeById(const std::string& id);

  bool addService(Service* a);
  const Service* getServiceById(const std::string& id) const;
  Service* getServiceById(const std::string& id);

 private:
  std::unordered_map<std::string, Agency*> _agencies;
  std::unordered_map<std::string, Stop*> _stops;
  std::unordered_map<std::string, Route*> _routes;
  std::unordered_map<std::string, Trip*> _trips;
  std::unordered_map<std::string, Shape*> _shapes;
  std::unordered_map<std::string, Service*> _services;
};
}}

#endif  // GTFSPARSER_GTFS_FEED_H_
