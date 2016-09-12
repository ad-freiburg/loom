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
typedef boost::geometry::model::polygon<Point> Polygon;
typedef boost::geometry::model::multi_polygon<Polygon> MultiPolygon;

// _____________________________________________________________________________
inline bool doubleEq(double a, double b) {
  return fabs(a - b) < 0.000001;
}

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
inline double angBetween(double p1x, double p1y, double q1x, double q1y) {
  double dY = q1y - p1y;
  double dX = q1x - p1x;
  return atan2(dY, dX) * (180.0 / M_PI);
}

// _____________________________________________________________________________
inline double angBetween(const Point& p1, const Point& q1) {
  return angBetween(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>());
}

// _____________________________________________________________________________
inline double dist(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1)*(x2 - x1)
              + (y2 - y1)*(y2 - y1));
}

// _____________________________________________________________________________
inline double innerProd(double x1, double y1, double x2, double y2,
    double x3, double y3) {
  double dx21 = x2-x1;
  double dx31 = x3-x1;
  double dy21 = y2-y1;
  double dy31 = y3-y1;
  double m12 = sqrt(dx21*dx21 + dy21*dy21);
  double m13 = sqrt(dx31*dx31 + dy31*dy31);
  double theta = acos((dx21*dx31 + dy21*dy31) / (m12 * m13));

  return theta * (180 / M_PI);
}

// _____________________________________________________________________________
inline double innerProd(const Point& a, const Point& b, const Point& c) {
  return innerProd(a.get<0>(), a.get<1>(), b.get<0>(), b.get<1>(), c.get<0>(), c.get<1>());
}

// _____________________________________________________________________________
inline double dist(const Point& p1, const Point& p2) {
  return dist(p1.get<0>(), p1.get<1>(), p2.get<0>(), p2.get<1>());
}

 // _____________________________________________________________________________
inline double distToSegment(double lax, double lay, double lbx, double lby,
                            double px, double py) {
  double d = dist(lax, lay, lbx, lby) * dist(lax, lay, lbx, lby);
  if (d == 0) return dist(px, py, lax, lay);

  double t = ((px - lax) * (lbx - lax)
    + (py - lay) * (lby - lay)) / d;

  if (t < 0) {
    return dist(px, py, lax, lay);
  } else if (t > 1) {
    return dist(px, py, lbx, lby);
  }

  return dist(px, py, lax + t * (lbx - lax), lay + t * (lby - lay));
}

// _____________________________________________________________________________
inline double distToSegment(const Point& la, const Point& lb, const Point& p) {
  return distToSegment(la.get<0>(), la.get<1>(), lb.get<0>(), lb.get<1>(), p.get<0>(),
    p.get<1>());
}

// _____________________________________________________________________________
inline Point projectOn(const Point& a, const Point& b, const Point& c) {
  if (doubleEq(a.get<0>(), b.get<0>()) && doubleEq(a.get<1>(), b.get<1>())) return a;
  if (doubleEq(a.get<0>(), c.get<0>()) && doubleEq(a.get<1>(), c.get<1>())) return a;
  if (doubleEq(b.get<0>(), c.get<0>()) && doubleEq(b.get<1>(), c.get<1>())) return b;

  double x;
  double y;

  if (c.get<0>() == a.get<0>()) {
    // infinite slope
    x = a.get<0>();
    y = b.get<1>();
  } else {
    double m = (double)(c.get<1>() - a.get<1>()) / (c.get<0>() - a.get<0>());
    double bb = (double)a.get<1>() - (m * a.get<0>());

    x = (m * b.get<1>() + b.get<0>() - m * bb) / (m * m + 1);
    y = (m * m * b.get<1>() + m * b.get<0>() + bb) / (m * m + 1);
  }

  Point ret = Point(x, y);

  bool isBetween = dist(a, c) > dist(a, ret) && dist(a, c) > dist(c, ret);
  bool nearer = dist(a, ret) < dist(c, ret);

  if (!isBetween) return nearer ? a : c;

  return ret;
}

}}}

#endif  // TRANSITMAP_UTIL_GEO_H_


