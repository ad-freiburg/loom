// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_CONFIG_GTFS2TOPOCONFIG_H_
#define GTFS2TOPO_CONFIG_GTFS2TOPOCONFIG_H_

#include <string>

namespace gtfs2topo {
namespace config {

struct Config {
  std::string projectionString;
  std::string inputFeedPath;

  bool ignoreGtfsDistances;
  bool ignoreDirections;

  size_t stationAggrLevel;

  double maxAggrDistance;
  double stationAggrDistance;

  uint8_t useMots;
};

}  // namespace config
}  // namespace gtfs2topo

#endif  // GTFS2TOPO_CONFIG_GTFS2TOPOCONFIG_H_
