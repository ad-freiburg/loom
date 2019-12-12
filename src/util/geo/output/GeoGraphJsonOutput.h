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

namespace util {
namespace geo {
namespace output {

class GeoGraphJsonOutput {
 public:
  inline GeoGraphJsonOutput(){};

  // print a graph to the provided path
  template <typename N, typename E>
  void print(const util::graph::Graph<N, E>& outG, std::ostream& str);

  // print a graph to the provided path, but treat coordinates as Web Mercator coordinates and reproject to WGS84
  template <typename N, typename E>
  void printLatLng(const util::graph::Graph<N, E>& outG, std::ostream& str);

 private:
  template <typename T>
  Line<T> createLine(const util::geo::Point<T>& a,
                     const util::geo::Point<T>& b);

  // print a graph to the provided path
  template <typename N, typename E>
  void printImpl(const util::graph::Graph<N, E>& outG, std::ostream& str, bool proj);
};

#include "util/geo/output/GeoGraphJsonOutput.tpp"
}
}
}

#endif  // UTIL_GEO_OUTPUT_GEOGRAPHJSONOUTPUT_H_
