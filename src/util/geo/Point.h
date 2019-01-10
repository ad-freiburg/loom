// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_POINT_H_
#define UTIL_GEO_POINT_H_

#include <vector>

namespace util {
namespace geo {

template <typename T>
class Point {
 public:
  Point() : _x(0), _y(0) {}
  Point(T x, T y) : _x(x), _y(y) {}
  T getX() const { return _x; }
  T getY() const { return _y; }

  void setX(T x) { _x = x; }
  void setY(T y) { _y = y; }

  Point<T> operator+(const Point<T>& p) const {
    return Point<T>(_x + p.getX(), _y + p.getY());
  }

  Point<T> operator-(const Point<T>& p) const {
    return Point<T>(_x - p.getX(), _y - p.getY());
  }

  bool operator==(const Point<T>& p) const {
    return p.getX() == _x && p.getY() == _y;
  }

  bool operator!=(const Point<T>& p) const {
    return !(*this == p);
  }

 private:
  T _x, _y;
};

template <typename T>
using MultiPoint = std::vector<Point<T>>;

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_POINT_H_
