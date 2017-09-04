// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_
#define UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_

#include <ostream>
#include <string>
#include "gtfs2topo/config/GraphBuilderConfig.h"
#include "gtfs2topo/graph/BuildGraph.h"
#include "json/json.hpp"
#include "util/String.h"
#include "gtfs2topo/graph/NodePL.h"
#include "gtfs2topo/graph/EdgePL.h"

using json = nlohmann::json;
using util::toString;

namespace gtfs2topo {

class GeoJsonOutput {
 public:
  GeoJsonOutput();
  template <typename N, typename E>
  void print(const util::graph::Graph<N, E>& outG);

 private:
};

#include "util/geo/output/GeoJsonOutput.tpp"

}

#endif  // UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_
