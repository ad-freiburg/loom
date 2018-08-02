// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_BOX_H_
#define UTIL_GEO_BOX_H_

#include "./Point.h"

namespace util {
namespace geo {

template <typename T>
class Box {
 public:
  // maximum inverse box as default value of box
  Box()
      : _ll(std::numeric_limits<T>::max(), std::numeric_limits<T>::max()),
        _ur(std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()) {}
  Box(const Point<T>& ll, const Point<T>& ur) : _ll(ll), _ur(ur) {}
  const Point<T>& getLowerLeft() const { return _ll; }
  const Point<T>& getUpperRight() const { return _ur; }

  Point<T>& getLowerLeft() { return _ll; }
  Point<T>& getUpperRight() { return _ur; }

  void setLowerLeft(const Point<T>& ll) { _ll = ll; }
  void setUpperRight(const Point<T>& ur) { _ur = ur; }

  bool operator==(const Box<T>& b) const {
    return getLowerLeft() == b.getLowerLeft() &&
           getUpperRight() == b.getUpperRight();
  }

  bool operator!=(const Box<T>& p) const { return !(*this == p); }

 private:
  Point<T> _ll, _ur;
};

template <typename T>
class RotatedBox {
 public:
  RotatedBox() : _box(), _deg(0), _center() {}
  RotatedBox(const Box<T>& box)
      : _box(box),
        _deg(0),
        _center(Point<T>(
            (box.getUpperRight().getX() - box.getLowerLeft().getX()) / T(2),
            (box.getUpperRight().getY() - box.getLowerLeft().getY()) / T(2))) {}
  RotatedBox(const Point<T>& ll, const Point<T>& ur)
      : _box(ll, ur),
        _deg(0),
        _center(Point<T>((ur.getX() - ll.getX()) / T(2),
                         (ur.getY() - ll.getY()) / T(2))) {}
  RotatedBox(const Box<T>& box, double deg)
      : _box(box),
        _deg(deg),
        _center(Point<T>(
            (box.getUpperRight().getX() - box.getLowerLeft().getX()) / T(2),
            (box.getUpperRight().getY() - box.getLowerLeft().getY()) / T(2))) {}
  RotatedBox(const Point<T>& ll, const Point<T>& ur, double deg)
      : _box(ll, ur),
        _deg(deg),
        _center(Point<T>((ur.getX() - ll.getX()) / T(2),
                         (ur.getY() - ll.getY()) / T(2))) {}
  RotatedBox(const Box<T>& box, double deg, const Point<T>& center)
      : _box(box), _deg(deg), _center(center) {}
  RotatedBox(const Point<T>& ll, const Point<T>& ur, double deg,
             const Point<T>& center)
      : _box(ll, ur), _deg(deg), _center(center) {}

  const Box<T>& getBox() const { return _box; }
  Box<T>& getBox() { return _box; }

  double getDegree() const { return _deg; }
  const Point<T>& getCenter() const { return _center; }
  Point<T>& getCenter() { return _center; }

  void setDegree(double deg) { _deg = deg; }

 private:
  Box<T> _box;
  double _deg;
  Point<T> _center;
};

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_BOX_H_
