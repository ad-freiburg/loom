// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_CONFIG_TOPOCONFIG_H_
#define TOPO_CONFIG_TOPOCONFIG_H_

#include <string>

namespace topo {
namespace config {

struct TopoConfig {
  double maxAggrDistance = 50;
  double maxLengthDev = 500;
  double maxTurnRestrCheckDist = 50;
  double turnInferFullTurnPen = 0;
  double fullTurnAngle = 0;
  double segmentLength = 5;
  bool outputStats = false;
  bool noInferRestrs = false;
  bool writeComponents = false;
  bool randomColors = false;
  double connectedCompDist = 10000;
  double smooth = 0;
  std::string componentsPath = "";
};

}  // namespace config
}  // namespace topo

#endif  // TOPO_CONFIG_TOPOCONFIG_H_
