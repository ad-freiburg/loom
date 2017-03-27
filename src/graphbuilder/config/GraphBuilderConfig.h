// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GRAPHBUILDER_CONFIG_GRAPHBUILDERCONFIG_H_
#define GRAPHBUILDER_CONFIG_GRAPHBUILDERCONFIG_H_

#include <string>

namespace graphbuilder {
namespace config {

struct Config {
  std::string projectionString;
  std::string inputFeedPath;
};

}  // namespace config
}  // namespace graphbuilder

#endif  // GRAPHBUILDER_CONFIG_GRAPHBUILDERCONFIG_H_
