// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GRAPHBUILDER_OUTPUT_JSONOUTPUT_H_
#define GRAPHBUILDER_OUTPUT_JSONOUTPUT_H_

#include <ostream>
#include <string>
#include "./../config/GraphBuilderConfig.h"
#include "./../graph/Graph.h"
#include "json/json.hpp"

using json = nlohmann::json;

namespace graphbuilder {

class JsonOutput {
 public:
  JsonOutput(const config::Config* cfg);
  void print(const graph::Graph& outG);

 private:
  const config::Config* _cfg;
};
}

#endif  // GRAPHBUILDER_OUTPUT_JSONOUTPUT_H_
