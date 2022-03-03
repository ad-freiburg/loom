// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef gtfs2graph_CONFIG_CONFIGREADER_H_
#define gtfs2graph_CONFIG_CONFIGREADER_H_

#include "gtfs2graph/config/GraphBuilderConfig.h"

namespace gtfs2graph {
namespace config {

class ConfigReader {
 public:
  ConfigReader();
  void read(Config* targetConfig, int argc, char** argv) const;
 private:
  void help(const char* bin) const;
};
}
}
#endif  // gtfs2graph_CONFIG_CONFIGREADER_H_
