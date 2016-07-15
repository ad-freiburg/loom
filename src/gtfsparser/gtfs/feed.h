// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_FEED_H_
#define GTFSPARSER_GTFS_FEED_H_

#include <vector>
#include <map>
#include <string>
#include "agency.h"
#include "stop.h"
#include "route.h"

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


 private:
  std::map<std::string, Agency*> _agencies;
  std::map<std::string, Stop*> _stops;
  std::map<std::string, Route*> _routes;
};
}}

#endif  // GTFSPARSER_GTFS_FEED_H_
