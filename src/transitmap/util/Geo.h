// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_UTIL_GEO_H_
#define TRANSITMAP_UTIL_GEO_H_

#define _USE_MATH_DEFINES

#include <boost/geometry.hpp>
#include <math.h>

// -------------------
// Geometry stuff
// ------------------

namespace transitmapper {
namespace util {
namespace geo {

// 2D cartesian coordinate
typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> Point;
typedef boost::geometry::model::linestring<Point> Line;
typedef boost::geometry::model::multi_linestring<Line> MultiLine;


// _____________________________________________________________________________
inline bool _onSegment(double px, double py,
                              double qx, double qy,
                              double rx, double ry) {
    if (qx <= std::max(px, rx) && qx >= std::min(px, rx) &&
        qy <= std::max(py, ry) && qy >= std::min(py, ry)) return true;

    return false;
}

// _____________________________________________________________________________
inline uint64_t _orientation(double px, double py,
                              double qx, double qy,
                              double rx, double ry) {
  int64_t val = (qy - py) * (rx - qx) - (qx - px) * (ry - qy);

  if (val == 0) return 0;
  return (val > 0) ? 1 : 2;
}

// _____________________________________________________________________________
inline bool intersects(double p1x, double p1y, double q1x, double q1y,
                        double p2x, double p2y, double q2x, double q2y) {
  int64_t o1 = _orientation(p1x, p1y, q1x, q1y, p2x, p2y);
  int64_t o2 = _orientation(p1x, p1y, q1x, q1y, q2x, q2y);
  int64_t o3 = _orientation(p2x, p2y, q2x, q2y, p1x, p1y);
  int64_t o4 = _orientation(p2x, p2y, q2x, q2y, q1x, q1y);

  return (o1 != o2 && o3 != o4)
    || (o1 == 0 && _onSegment(p1x, p1y, p2x, p2y, q1x, q1y))
    || (o2 == 0 && _onSegment(p1x, p1y, q2x, q2y, q1x, q1y))
    || (o3 == 0 && _onSegment(p2x, p2y, p1x, p1y, q2x, q2y))
    || (o4 == 0 && _onSegment(p2x, p2y, q1x, q1y, q2x, q2y));
}

// _____________________________________________________________________________
inline bool intersects(const Point& p1, const Point& q1, const Point& p2,
                        const Point& q2) {
  return intersects(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>(),
    p2.get<0>(), p2.get<1>(), q2.get<0>(), q2.get<1>());
}


// _____________________________________________________________________________
inline Point intersection(double p1x, double p1y, double q1x, double q1y,
                        double p2x, double p2y, double q2x, double q2y) {
  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));
  double u = (((q2x - p2x)*(p1y - p2y)) - ((q2y - p2y) * (p1x - p2x))) / a;

  return Point(p1x + (q1x - p1x)*u, p1y + (q1y - p1y) * u);
}

// _____________________________________________________________________________
inline Point intersection(const Point& p1, const Point& q1, const Point& p2,
                        const Point& q2) {
  return intersection(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>(),
    p2.get<0>(), p2.get<1>(), q2.get<0>(), q2.get<1>());
}

// _____________________________________________________________________________
inline bool lineIntersects(double p1x, double p1y, double q1x, double q1y,
                        double p2x, double p2y, double q2x, double q2y) {
  double EPSILON = 0.0000001;
  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));

  return a > EPSILON || a < -EPSILON;
}

// _____________________________________________________________________________
inline bool lineIntersects(const Point& p1, const Point& q1, const Point& p2,
                        const Point& q2) {
  return lineIntersects(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>(),
    p2.get<0>(), p2.get<1>(), q2.get<0>(), q2.get<1>());
}

// _____________________________________________________________________________
inline double segmentAngle(double p1x, double p1y, double q1x, double q1y) {
  double dY = q1y - p1y;
  double dX = q1x - p1x;
  return atan2(dY, dX) * (180.0 / M_PI);
}

// _____________________________________________________________________________
inline double segmentAngle(const Point& p1, const Point& q1) {
  return segmentAngle(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>());
}


}}}

#endif  // TRANSITMAP_UTIL_GEO_H_


