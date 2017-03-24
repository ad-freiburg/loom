// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PBUTIL_GEO_BEZIERCURVE_H_
#define PBUTIL_GEO_BEZIERCURVE_H_

#include <vector>
#include "./Geo.h"
#include "./PolyLine.h"
#include "./CubicSpline.h"

namespace pbutil {
namespace geo {

/**
 * Cubic spline class, not finished yet
 * (maybe not even the right choice and not needed...)
 */
class BezierCurve {

 public:
  BezierCurve(const Point& a, const Point& b, const Point& c,
    const Point& d);

  const PolyLine& render(double d);

 private:
  double _d;

  // the x and y polynoms for this spline
  CubicPolynom _xp;
  CubicPolynom _yp;

  // store the rendered polyline for quicker access
  PolyLine _rendered;
  bool _didRender;

  void recalcPolynoms(const Point& x, const Point& b,
    const Point& c, const Point& d);

  Point valueAt(double t) const;
};

}}

#endif  // PBUTIL_GEO_BEZIERCURVE_H_
