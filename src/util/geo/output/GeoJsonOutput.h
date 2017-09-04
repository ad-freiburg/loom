// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_
#define UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_

#include <ostream>
#include <string>
#include "json/json.hpp"
#include "util/String.h"
#include "util/graph/Graph.h"

using json = nlohmann::json;
using util::toString;
using util::graph::Graph;

namespace util {
namespace geo {
namespace output {

class GeoJsonOutput {
 public:
  GeoJsonOutput();
  template <typename N, typename E>
  void print(const Graph<N, E>& outG);

 private:
};

#include "util/geo/output/GeoJsonOutput.tpp"

}
}
}

#endif  // UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_
