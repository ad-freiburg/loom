// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename T>
CircularSegment<T>::CircularSegment(const Point<T>& a, double ang,
                            const Point<T>& c) : _a(a), _c(c), _renderD(0), _ang(ang)
     {

  _rad = dist(a, c);
  _s = fabs(_ang * _rad);
  _initAng = angBetween(c, a);
}

// _____________________________________________________________________________
template <typename T>
Point<T> CircularSegment<T>::valueAt(double ang) const {
  double xPos = _c.getX() + _rad * cos(ang);
  double yPos = _c.getY() + _rad * sin(ang);
  return Point<T>(xPos, yPos);
}

// _____________________________________________________________________________
template <typename T>
const PolyLine<T>& CircularSegment<T>::render(double d) {
  assert(d > 0);
  if (fabs(d - _renderD) < 0.001) return _rendered;
  _renderD = d;

  if (_s == 0) {
    _rendered << _a << _a;
    return _rendered;
  }

  _rendered.empty();
  double n = _s / d, dt = 1 / n, t = 0;

  bool cancel = false;
  while (true) {
    _rendered << valueAt(_initAng + t * _ang);
    t += dt;
    if (cancel) break;
    if (t > 1) {
      t = 1;
      cancel = true;
    }
  }

  return _rendered;
}
