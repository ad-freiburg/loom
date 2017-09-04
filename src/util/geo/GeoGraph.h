// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEOGRAPH_H_
#define UTIL_GEOGRAPH_H_

#include "util/geo/Geo.h"
#include "json/json.hpp"

using json = nlohmann::json;

namespace util {
namespace geograph {

class GeoEdgePL {
 public:
  virtual const util::geo::Line* getGeom() const = 0;
  virtual void getAttrs(json::object_t& obj) const = 0;
};

class GeoNodePL {
 public:
  virtual const util::geo::Point* getGeom() const = 0;
  virtual void getAttrs(json::object_t& obj) const = 0;
};

}  // namespace geograph
}  // namespace util

#endif  // UTIL_GEOGRAPH_H_
