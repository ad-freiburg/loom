// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_FEED_H_
#define GTFSPARSER_GTFS_FEED_H_

#include <unordered_map>
#include <string>
#include <iterator>
#include "Agency.h"
#include "Stop.h"
#include "Route.h"
#include "Trip.h"
#include "Shape.h"
#include "Service.h"

namespace gtfsparser {
namespace gtfs {

typedef std::unordered_map<std::string, Agency*> Agencies;
typedef std::unordered_map<std::string, Stop*> Stops;
typedef std::unordered_map<std::string, Route*> Routes;
typedef std::unordered_map<std::string, Trip*> Trips;
typedef std::unordered_map<std::string, Shape*> Shapes;
typedef std::unordered_map<std::string, Service*> Services;



class Feed {
 public:
  Feed() {}

  bool addAgency(Agency* a);
  const Agency* getAgencyById(const std::string& id) const;
  Agency* getAgencyById(const std::string& id);

  bool addStop(Stop* a);
  const Stop* getStopById(const std::string& id) const;
  Stop* getStopById(const std::string& id);
  Stops::const_iterator stopsBegin() const;
  Stops::const_iterator stopsEnd() const;
  Stops::iterator stopsBegin();
  Stops::iterator stopsEnd();

  bool addRoute(Route* a);
  const Route* getRouteById(const std::string& id) const;
  Route* getRouteById(const std::string& id);

  bool addTrip(Trip* a);
  const Trip* getTripById(const std::string& id) const;
  Trip* getTripById(const std::string& id);
  Trips::const_iterator tripsBegin() const;
  Trips::const_iterator tripsEnd() const;
  Trips::iterator tripsBegin();
  Trips::iterator tripsEnd();

  bool addShape(Shape* a);
  const Shape* getShapeById(const std::string& id) const;
  Shape* getShapeById(const std::string& id);
  Shapes::const_iterator shapesBegin() const;
  Shapes::const_iterator shapesEnd() const;
  Shapes::iterator shapesBegin();
  Shapes::iterator shapesEnd();

  bool addService(Service* a);
  const Service* getServiceById(const std::string& id) const;
  Service* getServiceById(const std::string& id);

 private:
  Agencies _agencies;
  Stops _stops;
  Routes _routes;
  Trips _trips;
  Shapes  _shapes;
  Services _services;
};
}  // namespace gtfs
}  // namespace gtfsparser

#endif  // GTFSPARSER_GTFS_FEED_H_
