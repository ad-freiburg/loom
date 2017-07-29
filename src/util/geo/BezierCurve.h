// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_BEZIERCURVE_H_
#define UTIL_GEO_BEZIERCURVE_H_

#include <vector>
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"

namespace util {
namespace geo {

struct CubicPolynom {
  CubicPolynom(double a, double b, double c, double d, double x)
      : a(a), b(b), c(c), d(d), x(x) {}
  CubicPolynom() : a(0), b(0), c(0), d(0), x(0) {}
  double a, b, c, d, x;

  double valueAt(double x) const;
};

/**
 * Bezier curve
 */
class BezierCurve {
 public:
  BezierCurve(const Point& a, const Point& b, const Point& c, const Point& d);

  const PolyLine& render(double d);

 private:
  double _d;

  // the x and y polynoms for this spline
  CubicPolynom _xp, _yp;

  // store the rendered polyline for quicker access
  PolyLine _rendered;
  bool _didRender;

  void recalcPolynoms(const Point& x, const Point& b, const Point& c,
                      const Point& d);

  Point valueAt(double t) const;
};
}
}

#endif  // UTIL_GEO_BEZIERCURVE_H_
