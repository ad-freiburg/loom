// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_CONFIG_OCTICONFIG_H_
#define OCTI_CONFIG_OCTICONFIG_H_

#include <string>
#include "octi/basegraph/BaseGraph.h"
#include "octi/basegraph/GridGraph.h"
#include "util/geo/Geo.h"

namespace octi {
namespace config {

enum OrderMethod {
  NUM_LINES = 0,
  LENGTH = 1,
  ADJ_ND_DEGREE = 2,
  ADJ_ND_LDEGREE = 3,
  GROWTH_DEG = 4,
  GROWTH_LDEG = 5,
  ALL = 99
};

struct Config {
  std::string gridSize = "100%";
  double borderRad = 45;

  std::string printMode = "linegraph";
  std::string optMode = "heur";
  std::string ilpPath;
  bool fromDot = false;
  bool deg2Heur = true;
  bool restrLocSearch = false;
  double enfGeoPen = 0;
  bool ilpNoSolve = false;
  int ilpTimeLimit = 60;
  int ilpNumThreads = 0;
  double ilpCacheThreshold = DBL_MAX;
  std::string ilpSolver = "gurobi";
  std::string ilpCacheDir = ".";

  bool skipOnError = false;
  bool retryOnError = false;

  double maxGrDist = 3;

  int heurLocSearchIters = 100;

  size_t abortAfter = -1;

  size_t hananIters = 1;
  bool writeStats = false;

  OrderMethod orderMethod;

  std::string obstaclePath;
  std::vector<util::geo::DPolygon> obstacles;

  octi::basegraph::BaseGraphType baseGraphType;

  octi::basegraph::Penalties pens;
};

}  // namespace config
}  // namespace octi

#endif  // OCTI_CONFIG_OCTICONFIG_H_
