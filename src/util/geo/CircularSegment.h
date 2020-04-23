// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_CIRCULARSEGMENT_H_
#define UTIL_GEO_CIRCULARSEGMENT_H_

#include <vector>
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"

namespace util {
namespace geo {

/**
 * Circular segment
 */
template <typename T>
class CircularSegment {
 public:
  CircularSegment(const Point<T>& a, double ang, const Point<T>& c);

  const PolyLine<T>& render(double d);

 private:
  // store the rendered polyline for quicker access
  PolyLine<T> _rendered;

  const Point<T>& _a, _c;
  double _renderD;

  double _ang, _rad, _s, _initAng;

  Point<T> valueAt(double t) const;
};

#include "util/geo/CircularSegment.tpp"
}
}

#endif  // UTIL_GEO_BEZIERCURVE_H_
