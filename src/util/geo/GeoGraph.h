// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEOGRAPH_H_
#define UTIL_GEOGRAPH_H_

#include <map>
#include "util/geo/Geo.h"
#include "util/json/Writer.h"

namespace util {
namespace geograph {

template<typename T>
class GeoEdgePL {
 public:
  virtual const util::geo::Line<T>* getGeom() const = 0;
  virtual json::Dict getAttrs() const = 0;
};

template<typename T>
class GeoNodePL {
 public:
  virtual const util::geo::Point<T>* getGeom() const = 0;
  virtual json::Dict getAttrs() const = 0;
};

}  // namespace geograph
}  // namespace util

#endif  // UTIL_GEOGRAPH_H_
