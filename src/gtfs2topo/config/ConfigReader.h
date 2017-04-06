// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_CONFIG_CONFIGREADER_H_
#define GTFS2TOPO_CONFIG_CONFIGREADER_H_

#include <boost/program_options.hpp>
#include <vector>
#include "./GraphBuilderConfig.h"

namespace gtfs2topo {
namespace config {

class ConfigReader {
 public:
  ConfigReader();
  void read(Config* targetConfig, int argc, char** argv) const;
};
}
}
#endif  // GTFS2TOPO_CONFIG_CONFIGREADER_H_
