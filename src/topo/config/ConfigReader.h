// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_CONFIG_CONFIGREADER_H_
#define TOPO_CONFIG_CONFIGREADER_H_

#include <boost/program_options.hpp>
#include <vector>
#include "topo/config/TopoConfig.h"

namespace topo {
namespace config {

class ConfigReader {
 public:
  ConfigReader();
  void read(TopoConfig* targetConfig, int argc, char** argv) const;
};
}
}
#endif  // TOPO_CONFIG_CONFIGREADER_H_
