// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_CONFIG_CONFIGREADER_H_
#define TRANSITMAP_CONFIG_CONFIGREADER_H_

#include <vector>
#include "transitmap/config/TransitMapConfig.h"

namespace transitmapper {
namespace config {

class ConfigReader {
 public:
  ConfigReader();
  void read(Config* targetConfig, int argc, char** argv) const;

 private:
  void help(const char* bin) const;
};
}  // namespace config
}  // namespace transitmapper
#endif  // TRANSITMAP_CONFIG_CONFIGREADER_H_
