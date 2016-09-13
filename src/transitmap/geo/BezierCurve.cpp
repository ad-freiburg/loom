// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include "./BezierCurve.h"
#include "./../util/Geo.h"
#include "./PolyLine.h"
#include "./CubicSpline.h"

using namespace transitmapper;
using namespace geo;

// _____________________________________________________________________________
BezierCurve::BezierCurve(const util::geo::Point& a, const util::geo::Point& b,
  const util::geo::Point& c, const util::geo::Point& d)
: _d(util::geo::dist(a, b)) {
  recalcPolynoms(a, b, c, d);
}

// _____________________________________________________________________________
void BezierCurve::recalcPolynoms(const util::geo::Point& a, const util::geo::Point& b,
  const util::geo::Point& c, const util::geo::Point& d) {

  _xp.a = a.get<0>();
  _xp.b = 3.0 * (b.get<0>() - a.get<0>());
  _xp.c = 3.0 * (c.get<0>() - b.get<0>()) - _xp.b;
  _xp.d = d.get<0>() - a.get<0>() - _xp.c - _xp.b;

  _yp.a = a.get<1>();
  _yp.b = 3.0 * (b.get<1>() - a.get<1>());
  _yp.c = 3.0 * (c.get<1>() - b.get<1>()) - _yp.b;
  _yp.d = d.get<1>() - a.get<1>() - _yp.c - _yp.b;

  _didRender = false;
}

// _____________________________________________________________________________
util::geo::Point BezierCurve::valueAt(double t) const {
  return util::geo::Point(_xp.valueAt(t), _yp.valueAt(t));
}

// _____________________________________________________________________________
const PolyLine& BezierCurve::render(double d) {
  if (_didRender) return _rendered;

  _rendered.empty();
  double n = _d / d;
  double dt = 1 / n;
  double t = 0;

  std::cout << dt << std::endl;

  bool cancel = false;
  while (true) {
    _rendered << valueAt(t);
    t += dt;
    if (cancel) break;
    if (t > 1) {
      t = 1;
      cancel = true;
    }
  }

  _didRender = true;
  return _rendered;
}
