// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_OUTPUT_GEOJSONOUTPUT_H_
#define GTFS2TOPO_OUTPUT_GEOJSONOUTPUT_H_

#include <ostream>
#include <string>
#include "gtfs2topo/config/GraphBuilderConfig.h"
#include "gtfs2topo/graph/BuildGraph.h"
#include "json/json.hpp"

using json = nlohmann::json;

namespace gtfs2topo {

class GeoJsonOutput {
 public:
  GeoJsonOutput(const config::Config* cfg);
  void print(const graph::BuildGraph& outG);

 private:
  const config::Config* _cfg;
};
}

#endif  // GTFS2TOPO_OUTPUT_GEOJSONOUTPUT_H_
