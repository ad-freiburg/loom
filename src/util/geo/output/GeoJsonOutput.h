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
#include "util/json/JsonWriter.h"

namespace util {
namespace geo {
namespace output {

typedef std::map<std::string, std::string> Attrs;

class GeoJsonOutput {
 public:
  GeoJsonOutput(std::ostream& str);
  ~GeoJsonOutput();
  template <typename T>
  void print(const Point<T>& p, Attrs attrs);
  template <typename T>
  void print(const Line<T>& l, Attrs attrs);
  void flush();

 private:
  json::JsonWriter _wr;
};

#include "util/geo/output/GeoJsonOutput.tpp"

}
}
}

#endif  // UTIL_GEO_OUTPUT_GEOJSONOUTPUT_H_
