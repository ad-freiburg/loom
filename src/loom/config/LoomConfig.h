// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_CONFIG_LOOMCONFIG_H_
#define LOOM_CONFIG_LOOMCONFIG_H_

#include <string>

namespace loom {
namespace config {

struct Config {
  std::string name;
  std::string outputPath;
  std::string dbgPath;

  std::string optimMethod = "comb-no-ilp";
  std::string MPSOutputPath;

  size_t optimRuns = 1;

  bool outOptGraph = false;

  bool outputStats = false;
  bool writeStats = false;
  bool pruneGraph = true;

  bool untangleGraph = true;
  bool fromDot = false;

  int ilpTimeLimit = -1;
  int ilpNumThreads = 0;

  double crossPenMultiSameSeg = 4;
  double crossPenMultiDiffSeg = 1;
  double separationPenWeight = 3;
  double stationCrossWeightSameSeg = 12;
  double stationCrossWeightDiffSeg = 3;
  double stationSeparationWeight = 9;

  std::string worldFilePath;

  std::string ilpSolver;
};

}  // namespace config
}  // namespace loom

#endif  // LOOM_CONFIG_LOOMCONFIG_H_
