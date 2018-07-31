// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_
#define UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_

#include <ostream>
#include <string>
#include <map>
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/json/Writer.h"

namespace util {
namespace geo {
namespace output {

class GeoJsonOutput {
 public:
  GeoJsonOutput(std::ostream& str);
  ~GeoJsonOutput();
  template <typename T>
  void print(const Point<T>& p, json::Val attrs);
  template <typename T>
  void print(const Line<T>& l, json::Val attrs);
  void flush();

 private:
  json::Writer _wr;
};

#include "util/geo/output/GeoJsonOutput.tpp"

}
}
}

#endif  // UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_
