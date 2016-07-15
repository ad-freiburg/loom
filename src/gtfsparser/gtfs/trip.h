// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef GTFSPARSER_GTFS_TRIP_H_
#define GTFSPARSER_GTFS_TRIP_H_

#include <stdint.h>
#include <string>
#include <sstream>
#include "route.h"

using std::exception;
using std::string;


namespace gtfsparser {
namespace gtfs {

class Trip {

 public:
  Trip() {};

  Route* getRoute() const {
    return _route;
  }

  // TODO: implement setters


 private:
  Route* _route;
};

}}

#endif  // GTFSPARSER_GTFS_TRIP_H_
