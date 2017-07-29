// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/BezierCurve.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"

using namespace util;
using namespace geo;

// _____________________________________________________________________________
BezierCurve::BezierCurve(const Point& a, const Point& b, const Point& c,
                         const Point& d)
    : _d(dist(a, d)) {
  assert(_d > 0);
  recalcPolynoms(a, b, c, d);
}

// _____________________________________________________________________________
void BezierCurve::recalcPolynoms(const Point& a, const Point& b, const Point& c,
                                 const Point& d) {
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
Point BezierCurve::valueAt(double t) const {
  return Point(_xp.valueAt(t), _yp.valueAt(t));
}

// _____________________________________________________________________________
const PolyLine& BezierCurve::render(double d) {
  assert(d > 0);
  if (_didRender) return _rendered;

  if (_d == 0) {
    _rendered << Point(_xp.a, _yp.a) << Point(_xp.a, _yp.a);
    return _rendered;
  }

  _rendered.empty();
  double n = _d / d, dt = 1 / n, t = 0;

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

// _____________________________________________________________________________
double CubicPolynom::valueAt(double atx) const {
  double dx = atx - x;
  return a + b * dx + c * dx * dx + d * dx * dx * dx;
}
