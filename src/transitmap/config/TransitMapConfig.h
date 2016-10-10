// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
#define TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

#include <string>

namespace transitmapper {
namespace config {

struct Config {

  double lineWidth;
  double lineSpacing;

  std::string projectionString;
  std::string inputFeedPath;
  std::string outputPath;

  std::string renderMethod;

  size_t optimIterations;

  double outputResolution;

  bool renderStations;
  bool renderNodeFronts;
  bool renderStationNames;

};

}  // namespace config
}  // namespace transitmapper

#endif  // TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
