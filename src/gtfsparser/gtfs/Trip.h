// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_TRIP_H_
#define GTFSPARSER_GTFS_TRIP_H_

#include <stdint.h>
#include <set>
#include <string>
#include "Route.h"
#include "Service.h"
#include "Shape.h"
#include "StopTime.h"
#include "Stop.h"

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

typedef std::set<StopTime, StopTimeCompare> StopTimes;


class Trip {

 public:
  enum WC_BIKE_ACCESSIBLE: uint8_t {
    NO_INFORMATION = 0,
    AT_LEAST_ONE = 1,
    NOT_ACCESSIBLE = 2
  };

  enum DIRECTION: uint8_t {
    OUTBOUND = 0,
    INBOUND = 1,
    NOT_SET = 2
  };

  Trip() {};
  Trip(std::string id, Route* r, Service* s, std::string hs,
    std::string short_name, DIRECTION dir, std::string blockid, Shape* shp,
    WC_BIKE_ACCESSIBLE wc, WC_BIKE_ACCESSIBLE ba);

  const std::string& getId() const;
  Route* getRoute() const;
  Service* getService() const;
  const std::string& getHeadsign() const;
  const std::string& getShortname() const;
  DIRECTION getDirection() const;
  const std::string& getBlockId() const;
  Shape* getShape() const;
  WC_BIKE_ACCESSIBLE getWheelchairAccessibility() const;
  WC_BIKE_ACCESSIBLE getBikesAllowed() const;
  const StopTimes& getStopTimes() const;
  bool addStopTime(const StopTime& t);

 private:
  std::string _id;
  Route* _route;
  Service* _service;
  std::string _headsign;
  std::string _short_name;
  DIRECTION _dir;
  std::string _block_id;
  Shape* _shape;
  WC_BIKE_ACCESSIBLE _wc;
  WC_BIKE_ACCESSIBLE _ba;

  StopTimes _stoptimes;
};

}}

#endif  // GTFSPARSER_GTFS_TRIP_H_
