// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GEO_BEZIERCURVE_H_
#define TRANSITMAP_GEO_BEZIERCURVE_H_

#include <vector>
#include "./../util/Geo.h"
#include "./PolyLine.h"
#include "./CubicSpline.h"

namespace transitmapper {
namespace geo {

/**
 * Cubic spline class, not finished yet
 * (maybe not even the right choice and not needed...)
 */
class BezierCurve {

 public:
  BezierCurve(const util::geo::Point& a, const util::geo::Point& b, const util::geo::Point& c,
    const util::geo::Point& d);

  const PolyLine& render(double d);

 private:
  double _d;

  // the x and y polynoms for this spline
  CubicPolynom _xp;
  CubicPolynom _yp;

  // store the rendered polyline for quicker access
  PolyLine _rendered;
  bool _didRender;

  void recalcPolynoms(const util::geo::Point& x, const util::geo::Point& b,
    const util::geo::Point& c, const util::geo::Point& d);

  util::geo::Point valueAt(double t) const;
};

}}

#endif  // TRANSITMAP_GEO_BEZIERCURVE_H_
