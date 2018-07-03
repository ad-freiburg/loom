// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_OUTPUT_GEOGRAPHJSONOUTPUT_H_
#define UTIL_GEO_OUTPUT_GEOGRAPHJSONOUTPUT_H_

#include <ostream>
#include <string>
#include "util/String.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/graph/Graph.h"

using util::toString;
using util::graph::Graph;

namespace util {
namespace geo {
namespace output {

class GeoGraphJsonOutput {
 public:
  inline GeoGraphJsonOutput(){};
  template <typename N, typename E>
  void print(const Graph<N, E>& outG, std::ostream& str);

 private:
  template <typename T>
  Line<T> createLine(const util::geo::Point<T>& a,
                     const util::geo::Point<T>& b);
};

#include "util/geo/output/GeoGraphJsonOutput.tpp"
}
}
}

#endif  // UTIL_GEO_OUTPUT_GEOGRAPHJSONOUTPUT_H_
