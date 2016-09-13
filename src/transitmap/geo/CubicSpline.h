// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GEO_CUBICSPLINE_H_
#define TRANSITMAP_GEO_CUBICSPLINE_H_

#include <vector>
#include "./PolyLine.h"

namespace transitmapper {
namespace geo {

struct CubicPolynom {
  CubicPolynom(double a, double b, double c, double d, double x)
  : a(a), b(b), c(c), d(d), x(x) {}
  CubicPolynom()
  : a(0), b(0), c(0), d(0), x(0) {}
  double a, b, c, d, x;

  double valueAt(double x) const;
};

/**
 * Cubic spline class, not finished yet
 * (maybe not even the right choice and not needed...)
 */
class CubicSpline {

 public:
  CubicSpline(const PolyLine& knots);
  CubicSpline(const util::geo::Line& knots);

  PolyLine render(double d) const;

 private:
  // the set of polynoms for this spline
  std::vector<CubicPolynom> _set;

  // store the rendered polyline for quicker access
  PolyLine _rendered;
  mutable bool _didRender;

  void findPolynoms(const util::geo::Line& knots, std::vector<CubicPolynom>* ret) const;
};

}}

#endif  // TRANSITMAP_GEO_CUBICSPLINE_H_
