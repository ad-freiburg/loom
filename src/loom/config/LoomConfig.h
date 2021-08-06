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

  std::string optimMethod;
  std::string MPSOutputPath;

  size_t optimRuns;

  bool separationOpt;
  bool outOptGraph;

  bool outputStats;
  bool renderStats;
  bool collapseLinePartners;

  bool createCoreOptimGraph;
  bool untangleGraph;
  bool simpleRenderForTwoEdgeNodes;

  int ilpTimeLimit;

  double crossPenMultiSameSeg;
  double crossPenMultiDiffSeg;
  double separationPenWeight;
  double stationCrossWeightSameSeg;
  double stationCrossWeightDiffSeg;
  double stationSeparationWeight;

  std::string worldFilePath;

  std::string ilpSolver;
};

}  // namespace config
}  // namespace loom

#endif  // LOOM_CONFIG_LOOMCONFIG_H_
