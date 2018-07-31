// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename T>
BezierCurve<T>::BezierCurve(const Point<T>& a, const Point<T>& b,
                            const Point<T>& c, const Point<T>& d)
    : _d(dist(a, d)) {
  assert(_d > 0);
  recalcPolynoms(a, b, c, d);
}

// _____________________________________________________________________________
template <typename T>
void BezierCurve<T>::recalcPolynoms(const Point<T>& a, const Point<T>& b,
                                    const Point<T>& c, const Point<T>& d) {
  _xp.a = a.getX();
  _xp.b = 3.0 * (b.getX() - a.getX());
  _xp.c = 3.0 * (c.getX() - b.getX()) - _xp.b;
  _xp.d = d.getX() - a.getX() - _xp.c - _xp.b;

  _yp.a = a.getY();
  _yp.b = 3.0 * (b.getY() - a.getY());
  _yp.c = 3.0 * (c.getY() - b.getY()) - _yp.b;
  _yp.d = d.getY() - a.getY() - _yp.c - _yp.b;

  _didRender = false;
}

// _____________________________________________________________________________
template <typename T>
Point<T> BezierCurve<T>::valueAt(double t) const {
  return Point<T>(_xp.valueAt(t), _yp.valueAt(t));
}

// _____________________________________________________________________________
template <typename T>
const PolyLine<T>& BezierCurve<T>::render(double d) {
  assert(d > 0);
  if (_didRender) return _rendered;

  if (_d == 0) {
    _rendered << Point<T>(_xp.a, _yp.a) << Point<T>(_xp.a, _yp.a);
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
