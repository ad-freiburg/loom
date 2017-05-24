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
  std::string outputPath;

  std::string renderMethod;
  std::string optimMethod;
  std::string glpkMPSOutputPath;
  std::string glpkHOutputPath;
  std::string glpkSolutionOutputPath;

  bool noOptim;
  bool splittingOpt;

  double outputResolution;
  double inputSmoothing;
  double innerGeometryPrecision;

  double inStationCrossPenalty;

  double outlineWidth;
  std::string outlineColor;

  bool renderStations;
  bool renderNodeFronts;
  bool renderNodeCircles;
  bool renderStationNames;

  bool useGlpkFeasibilityPump;
  bool useGlpkProximSearch;
  int glpkPSTimeLimit;
  int glpkTimeLimit;
};

}  // namespace config
}  // namespace transitmapper

#endif  // TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
