// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_CONFIG_TOPOCONFIG_H_
#define TOPO_CONFIG_TOPOCONFIG_H_

#include <string>

namespace topo {
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
}  // namespace topo

#endif  // TOPO_CONFIG_TOPOCONFIG_H_
