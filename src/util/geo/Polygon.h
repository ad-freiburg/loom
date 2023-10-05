// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_POLYGON_H_
#define UTIL_GEO_POLYGON_H_

#include <vector>
#include "./Box.h"
#include "./Line.h"
#include "./Point.h"

namespace util {
namespace geo {

template <typename T>
class Polygon {
 public:
  Polygon() {}

  Polygon(const Line<T>& l) : _outer(l) {}
  Polygon(const Box<T>& b)
      : _outer({b.getLowerLeft(),
                Point<T>(b.getLowerLeft().getX(), b.getUpperRight().getY()),
                b.getUpperRight(),
                Point<T>(b.getUpperRight().getX(), b.getLowerLeft().getY()),
                b.getLowerLeft()}) {}

  const Line<T>& getOuter() const { return _outer; }
  Line<T>& getOuter() { return _outer; }

  const std::vector<Line<T>>& getInners() const { return _inners; }
  std::vector<Line<T>>& getInners() { return _inners; }

 private:
  Line<T> _outer;
  std::vector<Line<T>> _inners;
};

template <typename T>
using MultiPolygon = std::vector<Polygon<T>>;

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_LINE_H_
