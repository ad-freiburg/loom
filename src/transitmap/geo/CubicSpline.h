// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GEO_CUBICSPLINE_H_
#define TRANSITMAP_GEO_CUBICSPLINE_H_

#include <vector>
#include "./PolyLine.h"

namespace transitmapper {
namespace geo {

struct Polynom {
  Polynom(double a, double b, double c, double d, double x)
  : a(a), b(b), c(c), d(d), x(x) {}
  double a, b, c, d, x;
};

class CubicSpline {

 public:
  CubicSpline(const PolyLine& knots);
  CubicSpline(const util::geo::Line& knots);

  PolyLine render(double d) const;

 private:
  // the set of polynoms for this spline
  std::vector<Polynom> _set;

  // store the rendered polyline for quicker access
  PolyLine _rendered;
  mutable bool _didRender;

  void findPolynoms(const util::geo::Line& knots, std::vector<Polynom>* ret) const;
};

}}

#endif  // TRANSITMAP_GEO_CUBICSPLINE_H_
