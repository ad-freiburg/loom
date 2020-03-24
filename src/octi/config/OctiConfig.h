// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_CONFIG_OCTICONFIG_H_
#define OCTI_CONFIG_OCTICONFIG_H_

#include <string>
#include "octi/basegraph/GridGraph.h"
#include "util/geo/Geo.h"

namespace octi {
namespace config {

struct Config {
  std::string gridSize;
  double borderRad;

  std::string printMode;
  std::string optMode;
  std::string ilpPath;
  bool fromDot;
  bool deg2Heur;
  bool restrLocSearch;
  double enfGeoPen;
  bool ilpNoSolve;

  double maxGrDist;

  std::string obstaclePath;
  std::vector<util::geo::DPolygon> obstacles;

  octi::basegraph::Penalties pens;
};

}  // namespace config
}  // namespace octi

#endif  // OCTI_CONFIG_OCTICONFIG_H_
