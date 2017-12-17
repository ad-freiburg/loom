// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_CONFIG_OCTICONFIG_H_
#define OCTI_CONFIG_OCTICONFIG_H_

#include <string>
#include "octi/gridgraph/GridGraph.h"

namespace octi {
namespace config {

struct Config {
  double gridCellWidth;

  std::string printMode;
  bool fromDot;

  octi::gridgraph::Penalties pens;
};

}  // namespace config
}  // namespace gtfs2topo

#endif  // OCTI_CONFIG_OCTICONFIG_H_
