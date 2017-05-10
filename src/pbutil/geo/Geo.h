// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_UTIL_GEO_H_
#define TRANSITMAP_UTIL_GEO_H_

#define _USE_MATH_DEFINES

#include <math.h>
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>

// -------------------
// Geometry stuff
// ------------------

namespace bgeo = boost::geometry;

namespace pbutil {
namespace geo {

// 2D cartesian coordinate
typedef bgeo::model::point<double, 2, bgeo::cs::cartesian> Point;
typedef bgeo::model::linestring<Point> Line;
typedef bgeo::model::multi_linestring<Line> MultiLine;
typedef bgeo::model::polygon<Point> Polygon;
typedef bgeo::model::multi_polygon<Polygon> MultiPolygon;
typedef bgeo::model::box<Point> Box;

// _____________________________________________________________________________
template <typename Geometry>
inline Geometry rotate(const Geometry& geo, double deg, const Point& center) {
  Geometry ret;

  bgeo::strategy::transform::translate_transformer<double, 2, 2> translate(
      -center.get<0>(), -center.get<1>());
  bgeo::strategy::transform::rotate_transformer<bgeo::degree, double, 2, 2>
      rotate(deg);
  bgeo::strategy::transform::translate_transformer<double, 2, 2> translateBack(
      center.get<0>(), center.get<1>());

  bgeo::strategy::transform::ublas_transformer<double, 2, 2> translateRotate(
      prod(rotate.matrix(), translate.matrix()));
  bgeo::strategy::transform::ublas_transformer<double, 2, 2> all(
      prod(translateBack.matrix(), translateRotate.matrix()));

  bgeo::transform(geo, ret, all);

  return ret;
}

// _____________________________________________________________________________
template <typename Geometry>
inline Geometry rotate(const Geometry& geo, double deg) {
  Point center;
  bgeo::centroid(geo, center);
  return rotate(geo, deg, center);
}

// _____________________________________________________________________________
inline bool doubleEq(double a, double b) { return fabs(a - b) < 0.000001; }

// _____________________________________________________________________________
inline bool intersects(const Point& p1, const Point& q1, const Point& p2,
                       const Point& q2) {
  /*
   * checks whether two line segments intersect
   */
  Line a;
  a.push_back(p1);
  a.push_back(q1);
  Line b;
  b.push_back(p2);
  b.push_back(q2);

  return bgeo::intersects(a, b);
}

// _____________________________________________________________________________
inline bool intersects(double p1x, double p1y, double q1x, double q1y,
                       double p2x, double p2y, double q2x, double q2y) {
  /*
   * checks whether two line segments intersect
   */
  Point p1(p1x, p1y);
  Point q1(q1x, q1y);
  Point p2(p2x, p2y);
  Point q2(q2x, q2y);

  return intersects(p1, q1, p2, q2);
}

// _____________________________________________________________________________
inline bool contains(const Point& p1, const Point& q1, const Point& p2,
                     const Point& q2) {
  Line a;
  a.push_back(p1);
  a.push_back(q1);
  Line b;
  b.push_back(p2);
  b.push_back(q2);

  return bgeo::within(a, b);
}

// _____________________________________________________________________________
inline bool contains(double p1x, double p1y, double q1x, double q1y, double p2x,
                     double p2y, double q2x, double q2y) {
  Point p1(p1x, p1y);
  Point q1(q1x, q1y);
  Point p2(p2x, p2y);
  Point q2(q2x, q2y);

  return intersects(p1, q1, p2, q2);
}

// _____________________________________________________________________________
inline Point intersection(double p1x, double p1y, double q1x, double q1y,
                          double p2x, double p2y, double q2x, double q2y) {
  /*
   * calculates the intersection between two line segments
   */
  if (p1x == q1x && p1y == q1y)
    return Point(p1x, p1y);  // TODO: <-- intersecting with a point??
  if (p2x == q1x && p2y == q1y) return Point(p2x, p2y);
  if (p2x == q2x && p2y == q2y)
    return Point(p2x, p2y);  // TODO: <-- intersecting with a point??
  if (p1x == q2x && p1y == q2y) return Point(p1x, p1y);

  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));
  double u = (((q2x - p2x) * (p1y - p2y)) - ((q2y - p2y) * (p1x - p2x))) / a;

  return Point(p1x + (q1x - p1x) * u, p1y + (q1y - p1y) * u);
}

// _____________________________________________________________________________
inline Point intersection(const Point& p1, const Point& q1, const Point& p2,
                          const Point& q2) {
  /*
   * calculates the intersection between two line segments
   */
  return intersection(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>(),
                      p2.get<0>(), p2.get<1>(), q2.get<0>(), q2.get<1>());
}

// _____________________________________________________________________________
inline bool lineIntersects(double p1x, double p1y, double q1x, double q1y,
                           double p2x, double p2y, double q2x, double q2y) {
  /*
   * checks whether two lines intersect
   */
  double EPSILON = 0.0000001;
  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));

  return a > EPSILON || a < -EPSILON;
}

// _____________________________________________________________________________
inline bool lineIntersects(const Point& p1, const Point& q1, const Point& p2,
                           const Point& q2) {
  /*
   * checks whether two lines intersect
   */
  return lineIntersects(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>(),
                        p2.get<0>(), p2.get<1>(), q2.get<0>(), q2.get<1>());
}

// _____________________________________________________________________________
inline double angBetween(double p1x, double p1y, double q1x, double q1y) {
  double dY = q1y - p1y;
  double dX = q1x - p1x;
  return atan2(dY, dX);
}

// _____________________________________________________________________________
inline double angBetween(const Point& p1, const Point& q1) {
  return angBetween(p1.get<0>(), p1.get<1>(), q1.get<0>(), q1.get<1>());
}

// _____________________________________________________________________________
inline double dist(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

// _____________________________________________________________________________
inline double innerProd(double x1, double y1, double x2, double y2, double x3,
                        double y3) {
  double dx21 = x2 - x1;
  double dx31 = x3 - x1;
  double dy21 = y2 - y1;
  double dy31 = y3 - y1;
  double m12 = sqrt(dx21 * dx21 + dy21 * dy21);
  double m13 = sqrt(dx31 * dx31 + dy31 * dy31);
  double theta = acos((dx21 * dx31 + dy21 * dy31) / (m12 * m13));

  return theta * (180 / M_PI);
}

// _____________________________________________________________________________
inline double innerProd(const Point& a, const Point& b, const Point& c) {
  return innerProd(a.get<0>(), a.get<1>(), b.get<0>(), b.get<1>(), c.get<0>(),
                   c.get<1>());
}

// _____________________________________________________________________________
inline double dist(const Point& p1, const Point& p2) {
  return bgeo::distance(p1, p2);
}

// _____________________________________________________________________________
inline double distToSegment(double lax, double lay, double lbx, double lby,
                            double px, double py) {
  double d = dist(lax, lay, lbx, lby) * dist(lax, lay, lbx, lby);
  if (d == 0) return dist(px, py, lax, lay);

  double t = ((px - lax) * (lbx - lax) + (py - lay) * (lby - lay)) / d;

  if (t < 0) {
    return dist(px, py, lax, lay);
  } else if (t > 1) {
    return dist(px, py, lbx, lby);
  }

  return dist(px, py, lax + t * (lbx - lax), lay + t * (lby - lay));
}

// _____________________________________________________________________________
inline double distToSegment(const Point& la, const Point& lb, const Point& p) {
  return distToSegment(la.get<0>(), la.get<1>(), lb.get<0>(), lb.get<1>(),
                       p.get<0>(), p.get<1>());
}

// _____________________________________________________________________________
inline Point projectOn(const Point& a, const Point& b, const Point& c) {
  if (doubleEq(a.get<0>(), b.get<0>()) && doubleEq(a.get<1>(), b.get<1>()))
    return a;
  if (doubleEq(a.get<0>(), c.get<0>()) && doubleEq(a.get<1>(), c.get<1>()))
    return a;
  if (doubleEq(b.get<0>(), c.get<0>()) && doubleEq(b.get<1>(), c.get<1>()))
    return b;

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

// _____________________________________________________________________________
inline double parallelity(const Box& box, const Line& line) {
  double ret = M_PI;

  double a = angBetween(box.min_corner(), Point(box.min_corner().get<0>(),
                                                box.max_corner().get<1>()));
  double b = angBetween(box.min_corner(), Point(box.max_corner().get<0>(),
                                                box.min_corner().get<1>()));
  double c = angBetween(box.max_corner(), Point(box.min_corner().get<0>(),
                                                box.max_corner().get<1>()));
  double d = angBetween(box.max_corner(), Point(box.max_corner().get<0>(),
                                                box.min_corner().get<1>()));

  double e = angBetween(line.front(), line.back());

  double vals[] = {a, b, c, d};

  for (double ang : vals) {
    double v = fabs(ang - e);
    if (v > M_PI) v = 2 * M_PI - v;
    if (v > M_PI / 2) v = M_PI - v;
    if (v < ret) ret = v;
  }

  return 1 - (ret / (M_PI / 4));
}

// _____________________________________________________________________________
inline double parallelity(const Box& box, const MultiLine& multiline) {
  double ret = 0;
  for (const Line& l : multiline) {
    ret += parallelity(box, l);
  }

  return ret / multiline.size();
}

// _____________________________________________________________________________
inline Polygon getOrientedEnvelope(Polygon pol) {
  // TODO: implement this nicer, works for now, but inefficient
  // see
  // https://geidav.wordpress.com/tag/gift-wrapping/#fn-1057-FreemanShapira1975
  // for a nicer algorithm

  Point center;
  bgeo::centroid(pol, center);

  Box tmpBox;
  bgeo::envelope(pol, tmpBox);
  double rotateDeg = 0;

  // rotate in 5 deg steps
  for (int i = 1; i < 360; i += 1) {
    pol = rotate(pol, 1, center);
    Box e;
    bgeo::envelope(pol, e);
    if (bgeo::area(tmpBox) > bgeo::area(e)) {
      tmpBox = e;
      rotateDeg = i;
    }
  }

  Polygon hull;
  bgeo::convex_hull(tmpBox, hull);
  return rotate(hull, -rotateDeg, center);
}

// _____________________________________________________________________________
inline Polygon getOrientedEnvelope(Polygon pol) {
  // TODO: implement this nicer, works for now, but inefficient
  // see
  // https://geidav.wordpress.com/tag/gift-wrapping/#fn-1057-FreemanShapira1975
  // for a nicer algorithm

  Point center;
  bgeo::centroid(pol, center);

  Box tmpBox;
  bgeo::envelope(pol, tmpBox);
  double rotateDeg = 0;

  // rotate in 5 deg steps
  for (int i = 1; i < 360; i += 1) {
    pol = rotate(pol, 1, center);
    Box e;
    bgeo::envelope(pol, e);
    if (bgeo::area(tmpBox) > bgeo::area(e)) {
      tmpBox = e;
      rotateDeg = i;
    }
  }

  Polygon hull;
  bgeo::convex_hull(tmpBox, hull);
  return rotate(hull, -rotateDeg, center);
}

// _____________________________________________________________________________
inline Polygon getOrientedEnvelope(const MultiLine& ml) {
  // first, get hull of multilines
  Polygon hull;
  bgeo::convex_hull(ml, hull);

  // get oriented envelope for hull
  hull = getOrientedEnvelope(hull);

  double minDeg;
  double score;

  for (double i = 0; i <= 360; i += .5) {

  }

}

}
}

#endif  // TRANSITMAP_UTIL_GEO_H_
