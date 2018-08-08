// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_GEO_H_
#define UTIL_GEO_GEO_H_

#define _USE_MATH_DEFINES

#include <math.h>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <sstream>
#include "util/Misc.h"
#include "util/geo/Box.h"
#include "util/geo/Line.h"
#include "util/geo/Point.h"
#include "util/geo/Polygon.h"

// -------------------
// Geometry stuff
// ------------------

namespace util {
namespace geo {

// convenience aliases

typedef Point<double> DPoint;
typedef Point<float> FPoint;
typedef Point<int> IPoint;

typedef LineSegment<double> DLineSegment;
typedef LineSegment<float> FLineSegment;
typedef LineSegment<int> ILineSegment;

typedef Line<double> DLine;
typedef Line<float> FLine;
typedef Line<int> ILine;

typedef Box<double> DBox;
typedef Box<float> FBox;
typedef Box<int> IBox;

typedef Polygon<double> DPolygon;
typedef Polygon<float> FPolygon;
typedef Polygon<int> IPolygon;

const static double EPSILON = 0.00000000001;
const static double RAD = 0.017453292519943295;  // PI/180

// _____________________________________________________________________________
template <typename T>
inline Box<T> pad(const Box<T>& box, double padding) {
  return Box<T>(Point<T>(box.getLowerLeft().getX() - padding,
                         box.getLowerLeft().getY() - padding),
                Point<T>(box.getUpperRight().getX() + padding,
                         box.getUpperRight().getY() + padding));
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> centroid(const Point<T> p) {
  return p;
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> centroid(const LineSegment<T> ls) {
  return Point<T>((ls.first.getX() + ls.second.getX()) / T(2),
                  (ls.first.getY() + ls.second.getY()) / T(2));
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> centroid(const Line<T> ls) {
  double x = 0, y = 0;
  for (const auto& p : ls) {
    x += p.getX();
    y += p.getY();
  }
  return Point<T>(x / T(ls.size()), y / T(ls.size()));
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> centroid(const Polygon<T> ls) {
  return centroid(ls.getOuter());
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> centroid(const Box<T> box) {
  return centroid(LineSegment<T>(box.getLowerLeft(), box.getUpperRight()));
}

// _____________________________________________________________________________
template <typename T, template <typename> class Geometry>
inline Point<T> centroid(std::vector<Geometry<T>> multigeo) {
  Line<T> a;
  for (const auto& g : multigeo) a.push_back(centroid(g));
  return centroid(a);
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> rotate(const Point<T>& p, double deg) {
  UNUSED(deg);
  return p;
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> rotate(Point<T> p, double deg, const Point<T>& c) {
  deg *= -RAD;
  double si = sin(deg);
  double co = cos(deg);
  p = p - c;

  return Point<T>(p.getX() * co - p.getY() * si,
                  p.getX() * si + p.getY() * co) +
         c;
}

// _____________________________________________________________________________
template <typename T>
inline LineSegment<T> rotate(LineSegment<T> geo, double deg,
                             const Point<T>& c) {
  geo.first = rotate(geo.first, deg, c);
  geo.second = rotate(geo.second, deg, c);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
inline LineSegment<T> rotate(LineSegment<T> geo, double deg) {
  return (geo, deg, centroid(geo));
}

// _____________________________________________________________________________
template <typename T>
inline Line<T> rotate(Line<T> geo, double deg, const Point<T>& c) {
  for (size_t i = 0; i < geo.size(); i++) geo[i] = rotate(geo[i], deg, c);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> rotate(Polygon<T> geo, double deg, const Point<T>& c) {
  for (size_t i = 0; i < geo.getOuter().size(); i++)
    geo.getOuter()[i] = rotate(geo.getOuter()[i], deg, c);
  return geo;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline std::vector<Geometry<T>> rotate(std::vector<Geometry<T>> multigeo,
                                       double deg, const Point<T>& c) {
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = rotate(multigeo[i], deg, c);
  return multigeo;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline std::vector<Geometry<T>> rotate(std::vector<Geometry<T>> multigeo,
                                       double deg) {
  auto c = centroid(multigeo);
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = rotate(multigeo[i], deg, c);
  return multigeo;
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> move(const Point<T>& geo, T x, T y) {
  return Point<T>(geo.getX() + x, geo.getY() + y);
}

// _____________________________________________________________________________
template <typename T>
inline Line<T> move(Line<T> geo, T x, T y) {
  for (size_t i = 0; i < geo.size(); i++) geo[i] = move(geo[i], x, y);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
inline LineSegment<T> move(LineSegment<T> geo, T x, T y) {
  geo.first = move(geo.first, x, y);
  geo.second = move(geo.second, x, y);
  return geo;
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> move(Polygon<T> geo, T x, T y) {
  for (size_t i = 0; i < geo.getOuter().size(); i++)
    geo.getOuter()[i] = move(geo.getOuter()[i], x, y);
  return geo;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline std::vector<Geometry<T>> move(std::vector<Geometry<T>> multigeo, T x,
                                     T y) {
  for (size_t i = 0; i < multigeo.size(); i++)
    multigeo[i] = move(multigeo[i], x, y);
  return multigeo;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> minbox() {
  return Box<T>();
}

// _____________________________________________________________________________
template <typename T>
inline RotatedBox<T> shrink(const RotatedBox<T>& b, double d) {
  double xd = b.getBox().getUpperRight().getX() - b.getBox().getLowerLeft().getX();
  double yd = b.getBox().getUpperRight().getY() - b.getBox().getLowerLeft().getY();

  if (xd <= 2 * d) d = xd / 2 - 1;
  if (yd <= 2 * d) d = yd / 2 - 1;

  Box<T> r(
      Point<T>(b.getBox().getLowerLeft().getX() + d, b.getBox().getLowerLeft().getY() + d),
      Point<T>(b.getBox().getUpperRight().getX() - d, b.getBox().getUpperRight().getY() - d));

  return RotatedBox<T>(r, b.getDegree(), b.getCenter());
}

// _____________________________________________________________________________
inline bool doubleEq(double a, double b) { return fabs(a - b) < 0.000001; }

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Point<T>& p, const Box<T>& box) {
  return p.getX() >= box.getLowerLeft().getX() &&
         p.getX() <= box.getUpperRight().getX() &&
         p.getY() >= box.getLowerLeft().getY() &&
         p.getY() <= box.getUpperRight().getY();
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Line<T>& l, const Box<T>& box) {
  for (const auto& p : l)
    if (!contains(p, box)) return false;
  return true;
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const LineSegment<T>& l, const Box<T>& box) {
  return contains(l.first, box) && contains(l.second, box);
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Box<T>& b, const Box<T>& box) {
  return contains(b.getLowerLeft(), box) && contains(b.getUpperRight(), box);
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Point<T>& p, const LineSegment<T>& ls) {
  return fabs(crossProd(p, ls)) < EPSILON && contains(p, getBoundingBox(ls));
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Point<T>& p, const Line<T>& l) {
  for (size_t i = 1; i < l.size(); i++) {
    if (contains(p, LineSegment<T>(l[i - 1], l[i]))) return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Point<T>& p, const Polygon<T>& poly) {
  // see https://de.wikipedia.org/wiki/Punkt-in-Polygon-Test_nach_Jordan
  int8_t c = -1;

  for (size_t i = 1; i < poly.getOuter().size(); i++) {
    c *= polyContCheck(p, poly.getOuter()[i - 1], poly.getOuter()[i]);
    if (c == 0) return true;
  }

  c *= polyContCheck(p, poly.getOuter().back(), poly.getOuter()[0]);

  return c >= 0;
}

// _____________________________________________________________________________
template <typename T>
inline int8_t polyContCheck(const Point<T>& a, Point<T> b, Point<T> c) {
  if (a.getY() == b.getY() && a.getY() == c.getY())
    return (!((b.getX() <= a.getX() && a.getX() <= c.getX()) ||
              (c.getX() <= a.getX() && a.getX() <= b.getX())));
  if (fabs(a.getY() - b.getY()) < EPSILON &&
      fabs(a.getX() - b.getX()) < EPSILON)
    return 0;
  if (b.getY() > c.getY()) {
    Point<T> tmp = b;
    b = c;
    c = tmp;
  }
  if (a.getY() <= b.getY() || a.getY() > c.getY()) {
    return 1;
  }

  double d = (b.getX() - a.getX()) * (c.getY() - a.getY()) -
             (b.getY() - a.getY()) * (c.getX() - a.getX());
  if (d > 0) return -1;
  if (d < 0) return 1;
  return 0;
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Polygon<T>& polyC, const Polygon<T>& poly) {
  for (const auto& p : polyC.getOuter()) {
    if (!contains(p, poly)) {
      return false;
    }
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Line<T>& l, const Polygon<T>& poly) {
  for (const auto& p : l) {
    if (!contains(p, poly)) {
      return false;
    }
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Line<T>& l, const Line<T>& other) {
  for (const auto& p : l) {
    if (!contains(p, other)) {
      return false;
    }
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Box<T>& b, const Polygon<T>& poly) {
  return contains(b.getLowerLeft(), poly) &&
         contains(b.getUpperRight(), poly) &&
         contains(Point<T>(b.getUpperRight().getX(), b.getLowerLeft().getY()),
                  poly) &&
         contains(Point<T>(b.getLowerLeft().getX(), b.getUpperRight().getY()),
                  poly);
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Polygon<T>& poly, const Box<T>& b) {
  for (const auto& p : poly.getOuter()) {
    if (!contains(p, b)) return false;
  }
  return true;
}

// _____________________________________________________________________________
template <typename T>
inline bool contains(const Polygon<T>& poly, const Line<T>& l) {
  for (const auto& p : poly.getOuter()) {
    if (!contains(p, l)) return false;
  }
  return true;
}

// _____________________________________________________________________________
template <template <typename> class GeometryA,
          template <typename> class GeometryB, typename T>
inline bool contains(const std::vector<GeometryA<T>>& multigeo,
                     const GeometryB<T>& geo) {
  for (const auto& g : multigeo)
    if (!contains(g, geo)) return false;
  return true;
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const LineSegment<T>& ls1, const LineSegment<T>& ls2) {
  return intersects(getBoundingBox(ls1), getBoundingBox(ls2)) &&
         (contains(ls1.first, ls2) || contains(ls1.second, ls2) ||
          contains(ls2.first, ls1) || contains(ls2.second, ls1) ||
          (((crossProd(ls1.first, ls2) < 0) ^
           (crossProd(ls1.second, ls2) < 0)) &&
          ((crossProd(ls2.first, ls1) < 0) ^ (crossProd(ls2.second, ls1) < 0))));
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Point<T>& a, const Point<T>& b, const Point<T>& c,
                       const Point<T>& d) {
  // legacy function
  return intersects(LineSegment<T>(a, b), LineSegment<T>(c, d));
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Line<T>& ls1, const Line<T>& ls2) {
  for (size_t i = 1; i < ls1.size(); i++) {
    for (size_t j = 1; j < ls2.size(); j++) {
      if (intersects(LineSegment<T>(ls1[i - 1], ls1[i]),
                     LineSegment<T>(ls2[j - 1], ls2[j])))
        return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Line<T>& l, const Point<T>& p) {
  return contains(l, p);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Point<T>& p, const Line<T>& l) {
  return intersects(l, p);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Polygon<T>& l, const Point<T>& p) {
  return contains(l, p);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Point<T>& p, const Polygon<T>& l) {
  return intersects(l, p);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Box<T>& b1, const Box<T>& b2) {
  return b1.getLowerLeft().getX() <= b2.getUpperRight().getX() &&
         b1.getUpperRight().getX() >= b2.getLowerLeft().getX() &&
         b1.getLowerLeft().getY() <= b2.getUpperRight().getY() &&
         b1.getUpperRight().getY() >= b2.getLowerLeft().getY();
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const LineSegment<T>& ls, const Box<T>& b) {
  if (intersects(ls, LineSegment<T>(b.getLowerLeft(),
                                    Point<T>(b.getUpperRight().getX(),
                                             b.getLowerLeft().getY()))))
    return true;
  if (intersects(ls, LineSegment<T>(b.getLowerLeft(),
                                    Point<T>(b.getLowerLeft().getX(),
                                             b.getUpperRight().getY()))))
    return true;
  if (intersects(ls, LineSegment<T>(b.getUpperRight(),
                                    Point<T>(b.getLowerLeft().getX(),
                                             b.getUpperRight().getY()))))
    return true;
  if (intersects(ls, LineSegment<T>(b.getUpperRight(),
                                    Point<T>(b.getUpperRight().getX(),
                                             b.getLowerLeft().getY()))))
    return true;

  return contains(ls, b);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const LineSegment<T>& ls, const Polygon<T>& p) {
  for (size_t i = 1; i < p.getOuter().size(); i++) {
    if (intersects(LineSegment<T>(p.getOuter()[i - 1], p.getOuter()[i]), ls))
      return true;
  }

  if (intersects(LineSegment<T>(p.getOuter().back(), p.getOuter().front()), ls))
    return true;

  return contains(ls, p);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Polygon<T>& p, const LineSegment<T>& ls) {
  return intersects(ls, p);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Box<T>& b, const LineSegment<T>& ls) {
  return intersects(ls, b);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Line<T>& l, const Box<T>& b) {
  for (size_t i = 1; i < l.size(); i++) {
    if (intersects(LineSegment<T>(l[i - 1], l[i]), b)) return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Box<T>& b, const Line<T>& l) {
  return intersects(l, b);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Point<T>& p, const Box<T>& b) {
  return contains(p, b);
}

// _____________________________________________________________________________
template <typename T>
inline bool intersects(const Box<T>& b, const Point<T>& p) {
  return intersects(p, b);
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> intersection(T p1x, T p1y, T q1x, T q1y, T p2x, T p2y, T q2x,
                             T q2y) {
  /*
   * calculates the intersection between two line segments
   */
  if (doubleEq(p1x, q1x) && doubleEq(p1y, q1y))
    return Point<T>(p1x, p1y);  // TODO: <-- intersecting with a point??
  if (doubleEq(p2x, q1x) && doubleEq(p2y, q1y)) return Point<T>(p2x, p2y);
  if (doubleEq(p2x, q2x) && doubleEq(p2y, q2y))
    return Point<T>(p2x, p2y);  // TODO: <-- intersecting with a point??
  if (doubleEq(p1x, q2x) && doubleEq(p1y, q2y)) return Point<T>(p1x, p1y);

  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));
  double u = (((q2x - p2x) * (p1y - p2y)) - ((q2y - p2y) * (p1x - p2x))) / a;

  return Point<T>(p1x + (q1x - p1x) * u, p1y + (q1y - p1y) * u);
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> intersection(const Point<T>& p1, const Point<T>& q1,
                             const Point<T>& p2, const Point<T>& q2) {
  /*
   * calculates the intersection between two line segments
   */
  return intersection(p1.getX(), p1.getY(), q1.getX(), q1.getY(), p2.getX(),
                      p2.getY(), q2.getX(), q2.getY());
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> intersection(const LineSegment<T>& s1,
                             const LineSegment<T>& s2) {
  return intersection(s1.first, s1.second, s2.first, s2.second);
}

// _____________________________________________________________________________
template <typename T>
inline bool lineIntersects(T p1x, T p1y, T q1x, T q1y, T p2x, T p2y, T q2x,
                           T q2y) {
  /*
   * checks whether two lines intersect
   */
  double a = ((q2y - p2y) * (q1x - p1x)) - ((q2x - p2x) * (q1y - p1y));

  return a > EPSILON || a < -EPSILON;
}

// _____________________________________________________________________________
template <typename T>
inline bool lineIntersects(const Point<T>& p1, const Point<T>& q1,
                           const Point<T>& p2, const Point<T>& q2) {
  /*
   * checks whether two lines intersect
   */
  return lineIntersects(p1.getX(), p1.getY(), q1.getX(), q1.getY(), p2.getX(),
                        p2.getY(), q2.getX(), q2.getY());
}

// _____________________________________________________________________________
inline double angBetween(double p1x, double p1y, double q1x, double q1y) {
  double dY = q1y - p1y;
  double dX = q1x - p1x;
  return atan2(dY, dX);
}

// _____________________________________________________________________________
template <typename T>
inline double angBetween(const Point<T>& p1, const Point<T>& q1) {
  return angBetween(p1.getX(), p1.getY(), q1.getX(), q1.getY());
}

// _____________________________________________________________________________
inline double dist(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const LineSegment<T>& ls, const Point<T>& p) {
  return distToSegment(ls, p);
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const Point<T>& p, const LineSegment<T>& ls) {
  return dist(ls, p);
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const LineSegment<T>& ls1, const LineSegment<T>& ls2) {
  if (intersects(ls1, ls2)) return 0;
  double d1 = dist(ls1.first, ls2);
  double d2 = dist(ls1.second, ls2);
  double d3 = dist(ls2.first, ls1);
  double d4 = dist(ls2.second, ls1);
  return std::min(d1, std::min(d2, (std::min(d3, d4))));
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const Point<T>& p, const Line<T>& l) {
  double d = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < l.size(); i++) {
    double dTmp = distToSegment(l[i-1], l[i], p);
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const Line<T>& l, const Point<T>& p) {
  return dist(p, l);
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const LineSegment<T>& ls, const Line<T>& l) {
  double d = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < l.size(); i++) {
    double dTmp = dist(ls, LineSegment<T>(l[i-1], l[i]));
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const Line<T>& l, const LineSegment<T>& ls) {
  return dist(ls, l);
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const Line<T>& la, const Line<T>& lb) {
  double d = std::numeric_limits<double>::infinity();
  for (size_t i = 1; i < la.size(); i++) {
    double dTmp = dist(LineSegment<T>(la[i-1], la[i]), lb);
    if (dTmp < EPSILON) return 0;
    if (dTmp < d) d = dTmp;
  }
  return d;
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
  double theta = acos(std::min((dx21 * dx31 + dy21 * dy31) / (m12 * m13), 1.0));

  return theta * (180 / M_PI);
}

// _____________________________________________________________________________
template <typename T>
inline double innerProd(const Point<T>& a, const Point<T>& b,
                        const Point<T>& c) {
  return innerProd(a.getX(), a.getY(), b.getX(), b.getY(), c.getX(), c.getY());
}

// _____________________________________________________________________________
template <typename T>
inline double crossProd(const Point<T>& a, const Point<T>& b) {
  return a.getX() * b.getY() - b.getX() * a.getY();
}

// _____________________________________________________________________________
template <typename T>
inline double crossProd(const Point<T>& p, const LineSegment<T>& ls) {
  return crossProd(
      Point<T>(ls.second.getX() - ls.first.getX(),
               ls.second.getY() - ls.first.getY()),
      Point<T>(p.getX() - ls.first.getX(), p.getY() - ls.first.getY()));
}

// _____________________________________________________________________________
template <typename T>
inline double dist(const Point<T>& p1, const Point<T>& p2) {
  return dist(p1.getX(), p1.getY(), p2.getX(), p2.getY());
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const Point<T>& p) {
  std::stringstream ss;
  ss << "POINT (" << p.getX() << " " << p.getY() << ")";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const std::vector<Point<T>>& p) {
  std::stringstream ss;
  ss << "MULTIPOINT (";
  for (size_t i = 0; i < p.size(); i++) {
    if (i) ss << ", ";
    ss << "(" << p[i].getX() << " " << p[i].getY() << ")";
  }
  ss << ")";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const Line<T>& l) {
  std::stringstream ss;
  ss << "LINESTRING (";
  for (size_t i = 0; i < l.size(); i++) {
    if (i) ss << ", ";
    ss << l[i].getX() << " " << l[i].getY();
  }
  ss << ")";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const std::vector<Line<T>>& ls) {
  std::stringstream ss;
  ss << "MULTILINESTRING (";

  for (size_t j = 0; j < ls.size(); j++) {
    if (j) ss << ", ";
    ss << "(";
    for (size_t i = 0; i < ls[j].size(); i++) {
      if (i) ss << ", ";
      ss << ls[j][i].getX() << " " << ls[j][i].getY();
    }
    ss << ")";
  }

  ss << ")";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const LineSegment<T>& l) {
  return getWKT(Line<T>{l.first, l.second});
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const Box<T>& l) {
  std::stringstream ss;
  ss << "POLYGON ((";
  ss << l.getLowerLeft().getX() << " " << l.getLowerLeft().getY();
  ss << ", " << l.getUpperRight().getX() << " " << l.getLowerLeft().getY();
  ss << ", " << l.getUpperRight().getX() << " " << l.getUpperRight().getY();
  ss << ", " << l.getLowerLeft().getX() << " " << l.getUpperRight().getY();
  ss << ", " << l.getLowerLeft().getX() << " " << l.getLowerLeft().getY();
  ss << "))";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const Polygon<T>& p) {
  std::stringstream ss;
  ss << "POLYGON ((";
  for (size_t i = 0; i < p.getOuter().size(); i++) {
    ss << p.getOuter()[i].getX() << " " << p.getOuter()[i].getY() << ", ";
  }
  ss << p.getOuter().front().getX() << " " << p.getOuter().front().getY();
  ss << "))";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
inline std::string getWKT(const std::vector<Polygon<T>>& ls) {
  std::stringstream ss;
  ss << "MULTIPOLYGON (";

  for (size_t j = 0; j < ls.size(); j++) {
    if (j) ss << ", ";
    ss << "((";
    for (size_t i = 0; i < ls[j].getOuter().size(); i++) {
      ss << ls[j].getOuter()[i].getX() << " " << ls[j].getOuter()[i].getY()
         << ", ";
    }
    ss << ls[j].getOuter().front().getX() << " "
       << ls[j].getOuter().front().getY();
    ss << "))";
  }

  ss << ")";
  return ss.str();
}

// _____________________________________________________________________________
template <typename T>
inline double len(const Point<T>& g) {
  UNUSED(g);
  return 0;
}

// _____________________________________________________________________________
template <typename T>
inline double len(const Line<T>& g) {
  double ret = 0;
  for (size_t i = 1; i < g.size(); i++) ret += dist(g[i - 1], g[i]);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> simplify(const Point<T>& g, double d) {
  UNUSED(d);
  return g;
}

// _____________________________________________________________________________
template <typename T>
inline LineSegment<T> simplify(const LineSegment<T>& g, double d) {
  UNUSED(d);
  return g;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> simplify(const Box<T>& g, double d) {
  UNUSED(d);
  return g;
}

// _____________________________________________________________________________
template <typename T>
inline RotatedBox<T> simplify(const RotatedBox<T>& g, double d) {
  UNUSED(d);
  return g;
}

// _____________________________________________________________________________
template <typename T>
inline Line<T> simplify(const Line<T>& g, double d) {
  // douglas peucker
  double maxd = 0;
  size_t maxi = 0;
  for (size_t i = 1; i < g.size() - 1; i++) {
    double dt = distToSegment(g.front(), g.back(), g[i]);
    if (dt > maxd) {
      maxi = i;
      maxd = dt;
    }
  }

  if (maxd > d) {
    auto a = simplify(Line<T>(g.begin(), g.begin() + maxi + 1), d);
    const auto& b = simplify(Line<T>(g.begin() + maxi, g.end()), d);
    a.insert(a.end(), b.begin() + 1, b.end());

    return a;
  }

  return Line<T>{g.front(), g.back()};
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> simplify(const Polygon<T>& g, double d) {
  auto simple = simplify(g, d);
  std::rotate(simple.begin(), simple.begin() + simple.size() / 2, simple.end());
  simple = simplify(simple, d);
  return Polygon<T>(simple);
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
template <typename T>
inline double distToSegment(const Point<T>& la, const Point<T>& lb,
                            const Point<T>& p) {
  return distToSegment(la.getX(), la.getY(), lb.getX(), lb.getY(), p.getX(),
                       p.getY());
}

// _____________________________________________________________________________
template <typename T>
inline double distToSegment(const LineSegment<T>& ls, const Point<T>& p) {
  return distToSegment(ls.first.getX(), ls.first.getY(), ls.second.getX(),
                       ls.second.getY(), p.getX(), p.getY());
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> projectOn(const Point<T>& a, const Point<T>& b,
                          const Point<T>& c) {
  if (doubleEq(a.getX(), b.getX()) && doubleEq(a.getY(), b.getY())) return a;
  if (doubleEq(a.getX(), c.getX()) && doubleEq(a.getY(), c.getY())) return a;
  if (doubleEq(b.getX(), c.getX()) && doubleEq(b.getY(), c.getY())) return b;

  T x, y;

  if (c.getX() == a.getX()) {
    // infinite slope
    x = a.getX();
    y = b.getY();
  } else {
    double m = (double)(c.getY() - a.getY()) / (c.getX() - a.getX());
    double bb = (double)a.getY() - (m * a.getX());

    x = (m * b.getY() + b.getX() - m * bb) / (m * m + 1.0);
    y = (m * m * b.getY() + m * b.getX() + bb) / (m * m + 1.0);
  }

  Point<T> ret = Point<T>(x, y);

  bool isBetween = dist(a, c) > dist(a, ret) && dist(a, c) > dist(c, ret);
  bool nearer = dist(a, ret) < dist(c, ret);

  if (!isBetween) return nearer ? a : c;

  return ret;
}

// _____________________________________________________________________________
template <typename T>
inline double parallelity(const Box<T>& box, const Line<T>& line) {
  double ret = M_PI;

  double a = angBetween(
      box.getLowerLeft(),
      Point<T>(box.getLowerLeft().getX(), box.getUpperRight().getY()));
  double b = angBetween(
      box.getLowerLeft(),
      Point<T>(box.getUpperRight().getX(), box.getLowerLeft().getY()));
  double c = angBetween(
      box.getUpperRight(),
      Point<T>(box.getLowerLeft().getX(), box.getUpperRight().getY()));
  double d = angBetween(
      box.getUpperRight(),
      Point<T>(box.getUpperRight().getX(), box.getLowerLeft().getY()));

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
template <typename T>
inline double parallelity(const Box<T>& box, const MultiLine<T>& multiline) {
  double ret = 0;
  for (const Line<T>& l : multiline) {
    ret += parallelity(box, l);
  }

  return ret / static_cast<float>(multiline.size());
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline RotatedBox<T> getOrientedEnvelope(Geometry<T> pol) {
  // TODO: implement this nicer, works for now, but inefficient
  // see
  // https://geidav.wordpress.com/tag/gift-wrapping/#fn-1057-FreemanShapira1975
  // for a nicer algorithm

  Point<T> center = centroid(pol);
  Box<T> tmpBox = getBoundingBox(pol);
  double rotateDeg = 0;

  // rotate in 1 deg steps
  for (int i = 1; i < 360; i += 1) {
    pol = rotate(pol, 1, center);
    Box<T> e = getBoundingBox(pol);
    if (area(tmpBox) > area(e)) {
      tmpBox = e;
      rotateDeg = i;
    }
  }

  return RotatedBox<T>(tmpBox, -rotateDeg, center);
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> extendBox(const Box<T>& a, Box<T> b) {
  b = extendBox(a.getLowerLeft(), b);
  b = extendBox(a.getUpperRight(), b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> extendBox(const Point<T>& p, Box<T> b) {
  if (p.getX() < b.getLowerLeft().getX()) b.getLowerLeft().setX(p.getX());
  if (p.getY() < b.getLowerLeft().getY()) b.getLowerLeft().setY(p.getY());

  if (p.getX() > b.getUpperRight().getX()) b.getUpperRight().setX(p.getX());
  if (p.getY() > b.getUpperRight().getY()) b.getUpperRight().setY(p.getY());
  return b;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> getBoundingBox(const Point<T>& p) {
  return Box<T>(p, p);
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> getBoundingBox(const Line<T>& l) {
  Box<T> ret;
  for (const auto& p : l) ret = extendBox(p, ret);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> getBoundingBox(const Polygon<T>& pol) {
  Box<T> ret;
  for (const auto& p : pol.getOuter()) ret = extendBox(p, ret);
  return ret;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> getBoundingBox(const LineSegment<T>& ls) {
  Box<T> b;
  b = extendBox(ls.first, b);
  b = extendBox(ls.second, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> getBoundingBox(const Box<T>& b) {
  return b;
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline Box<T> getBoundingBox(const std::vector<Geometry<T>>& multigeo) {
  Box<T> b;
  b = extendBox(multigeo, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> convexHull(const Point<T>& p) {
  return Polygon<T>({p});
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> convexHull(const Box<T>& b) {
  return Polygon<T>(b);
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> convexHull(const LineSegment<T>& b) {
  return Polygon<T>(Line<T>{b.first, b.second});
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> convexHull(const RotatedBox<T>& b) {
  auto p = convexHull(b.getBox());
  p = rotate(p, b.getDegree(), b.getCenter());
  return p;
}

// _____________________________________________________________________________
template <typename T>
inline size_t convexHullImpl(const MultiPoint<T>& a, size_t p1, size_t p2,
                             Line<T>* h, uint8_t d) {
  // quickhull by Barber, Dobkin & Huhdanpaa
  Point<T> pa;
  bool found = false;
  double maxDist = 0;
  for (const auto& p : a) {
    double tmpDist = distToSegment((*h)[p1], (*h)[p2], p);
    double cp = crossProd(p, LineSegment<T>((*h)[p1], (*h)[p2]));
    if (((cp > 0 && !d) || (cp < 0 && d)) && tmpDist > maxDist) {
      pa = p;
      found = true;
      maxDist = tmpDist;
    }
  }

  if (!found) return 0;

  h->insert(h->begin() + p2 + !d, pa);
  size_t in = 1 + convexHullImpl(a, p1, p2 + !d, h, d);
  return in + convexHullImpl(a, p2 + in * d + 1 - 2 * d, p2 + in * d, h, d);
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> convexHull(const MultiPoint<T>& l) {
  if (l.size() == 2) return convexHull(LineSegment<T>(l[0], l[1]));
  if (l.size() == 1) return convexHull(l[0]);

  Point<T> left(std::numeric_limits<T>::max(), 0);
  Point<T> right(std::numeric_limits<T>::lowest(), 0);
  for (const auto& p : l) {
    if (p.getX() <= left.getX()) left = p;
    if (p.getX() >= right.getX()) right = p;
  }

  Line<T> hull{left, right};
  convexHullImpl(l, 0, 1, &hull, 1);
  convexHullImpl(l, 0, hull.size() - 1, &hull, 0);
  return Polygon<T>(hull);
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> convexHull(const Polygon<T>& p) {
  return convexHull(p.getOuter());
}

// _____________________________________________________________________________
template <typename T>
inline Polygon<T> convexHull(const MultiLine<T>& ls) {
  MultiPoint<T> mp;
  for (const auto& l : ls) mp.insert(mp.end(), l.begin(), l.end());
  return convexHull(mp);
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> extendBox(const Line<T>& l, Box<T> b) {
  for (const auto& p : l) b = extendBox(p, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> extendBox(const LineSegment<T>& ls, Box<T> b) {
  b = extendBox(ls.first, b);
  b = extendBox(ls.second, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
inline Box<T> extendBox(const Polygon<T>& ls, Box<T> b) {
  return extendBox(ls.getOuter(), b);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline Box<T> extendBox(const std::vector<Geometry<T>>& multigeom, Box<T> b) {
  for (const auto& g : multigeom) b = extendBox(g, b);
  return b;
}

// _____________________________________________________________________________
template <typename T>
inline double area(const Point<T>& b) {
  UNUSED(b);
  return 0;
}

// _____________________________________________________________________________
template <typename T>
inline double area(const LineSegment<T>& b) {
  UNUSED(b);
  return 0;
}

// _____________________________________________________________________________
template <typename T>
inline double area(const Line<T>& b) {
  UNUSED(b);
  return 0;
}

// _____________________________________________________________________________
template <typename T>
inline double area(const Box<T>& b) {
  return (b.getUpperRight().getX() - b.getLowerLeft().getX()) *
         (b.getUpperRight().getY() - b.getLowerLeft().getY());
}

// _____________________________________________________________________________
template <typename T>
inline double area(const Polygon<T>& b) {
  double ret = 0;
  size_t j = b.getOuter().size() - 1;
  for (size_t i = 0; i < b.getOuter().size(); i++) {
    ret += (b.getOuter()[j].getX() + b.getOuter()[i].getX()) *
           (b.getOuter()[j].getY() - b.getOuter()[i].getY());
    j = i;
  }

  return fabs(ret / 2.0);
}

// _____________________________________________________________________________
template <typename T>
inline double commonArea(const Box<T>& ba, const Box<T>& bb) {
  T l = std::max(ba.getLowerLeft().getX(), bb.getLowerLeft().getX());
  T r = std::min(ba.getUpperRight().getX(), bb.getUpperRight().getX());
  T b = std::max(ba.getLowerLeft().getY(), bb.getLowerLeft().getY());
  T t = std::min(ba.getUpperRight().getY(), bb.getUpperRight().getY());

  if (l > r || b > t) return 0;
  return (r - l) * (t - b);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline RotatedBox<T> getFullEnvelope(std::vector<Geometry<T>> pol) {
  Point<T> center = centroid(pol);
  Box<T> tmpBox = getBoundingBox(pol);
  double rotateDeg = 0;

  std::vector<Polygon<T>> ml;

  // rotate in 5 deg steps
  for (int i = 1; i < 360; i += 1) {
    pol = rotate(pol, 1, center);
    Polygon<T> hull = convexHull(pol);
    ml.push_back(hull);
    Box<T> e = getBoundingBox(pol);
    if (area(tmpBox) > area(e)) {
      tmpBox = e;
      rotateDeg = i;
    }
  }

  tmpBox = getBoundingBox(ml);

  return RotatedBox<T>(tmpBox, rotateDeg, center);
}

// _____________________________________________________________________________
template <template <typename> class Geometry, typename T>
inline RotatedBox<T> getFullEnvelope(const Geometry<T> pol) {
  std::vector<Geometry<T>> mult;
  mult.push_back(pol);
  return getFullEnvelope(mult);
}

// _____________________________________________________________________________
template <typename T>
inline RotatedBox<T> getOrientedEnvelopeAvg(MultiLine<T> ml) {
  MultiLine<T> orig = ml;
  // get oriented envelope for hull
  RotatedBox<T> rbox = getFullEnvelope(ml);
  Point<T> center = centroid(rbox.getBox());

  ml = rotate(ml, -rbox.getDegree() - 45, center);

  double bestDeg = -45;
  double score = parallelity(rbox.getBox(), ml);

  for (double i = -45; i <= 45; i += .5) {
    ml = rotate(ml, -.5, center);
    double p = parallelity(rbox.getBox(), ml);
    if (parallelity(rbox.getBox(), ml) > score) {
      bestDeg = i;
      score = p;
    }
  }

  rbox.setDegree(rbox.getDegree() + bestDeg);

  // move the box along 45deg angles from its origin until it fits the ml
  // = until the intersection of its hull and the box is largest
  Polygon<T> p = convexHull(rbox);
  p = rotate(p, -rbox.getDegree(), rbox.getCenter());

  Polygon<T> hull = convexHull(orig);
  hull = rotate(hull, -rbox.getDegree(), rbox.getCenter());

  Box<T> box = getBoundingBox(hull);
  rbox = RotatedBox<T>(box, rbox.getDegree(), rbox.getCenter());

  return rbox;
}

// _____________________________________________________________________________
template <typename T>
inline Line<T> densify(const Line<T>& l, double d) {
  if (!l.size()) return l;

  Line<T> ret;
  ret.reserve(l.size());
  ret.push_back(l.front());

  for (size_t i = 1; i < l.size(); i++) {
    double segd = dist(l[i - 1], l[i]);
    double dx = (l[i].getX() - l[i - 1].getX()) / segd;
    double dy = (l[i].getY() - l[i - 1].getY()) / segd;
    double curd = d;
    while (curd < segd) {
      ret.push_back(
          Point<T>(l[i - 1].getX() + dx * curd, l[i - 1].getY() + dy * curd));
      curd += d;
    }

    ret.push_back(l[i]);
  }

  return ret;
}

// _____________________________________________________________________________
template <typename T>
inline double frechetDistC(size_t i, size_t j, const Line<T>& p,
                           const Line<T>& q,
                           std::vector<std::vector<double>>& ca) {
  // based on Eiter / Mannila
  // http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

  if (ca[i][j] > -1)
    return ca[i][j];
  else if (i == 0 && j == 0)
    ca[i][j] = dist(p[0], q[0]);
  else if (i > 0 && j == 0)
    ca[i][j] = std::max(frechetDistC(i - 1, 0, p, q, ca), dist(p[i], q[0]));
  else if (i == 0 && j > 0)
    ca[i][j] = std::max(frechetDistC(0, j - 1, p, q, ca), dist(p[0], q[j]));
  else if (i > 0 && j > 0)
    ca[i][j] = std::max(std::min(std::min(frechetDistC(i - 1, j, p, q, ca),
                                          frechetDistC(i - 1, j - 1, p, q, ca)),
                                 frechetDistC(i, j - 1, p, q, ca)),
                        dist(p[i], q[j]));
  else
    ca[i][j] = std::numeric_limits<double>::infinity();

  return ca[i][j];
}

// _____________________________________________________________________________
template <typename T>
inline double frechetDist(const Line<T>& a, const Line<T>& b, double d) {
  // based on Eiter / Mannila
  // http://www.kr.tuwien.ac.at/staff/eiter/et-archive/cdtr9464.pdf

  auto p = densify(a, d);
  auto q = densify(b, d);

  std::vector<std::vector<double>> ca(p.size(),
                                      std::vector<double>(q.size(), -1.0));
  double fd = frechetDistC(p.size() - 1, q.size() - 1, p, q, ca);

  return fd;
}

// _____________________________________________________________________________
template <typename T>
inline double accFrechetDistC(const Line<T>& a, const Line<T>& b, double d) {
  auto p = densify(a, d);
  auto q = densify(b, d);

  std::vector<std::vector<double>> ca(p.size(),
                                      std::vector<double>(q.size(), 0));

  for (size_t i = 0; i < p.size(); i++)
    ca[i][0] = std::numeric_limits<double>::infinity();
  for (size_t j = 0; j < q.size(); j++)
    ca[0][j] = std::numeric_limits<double>::infinity();
  ca[0][0] = 0;

  for (size_t i = 1; i < p.size(); i++) {
    for (size_t j = 1; j < q.size(); j++) {
      double d = util::geo::dist(p[i], q[j]) * util::geo::dist(p[i], p[i - 1]);
      ca[i][j] =
          d + std::min(ca[i - 1][j], std::min(ca[i][j - 1], ca[i - 1][j - 1]));
    }
  }

  return ca[p.size() - 1][q.size() - 1];
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> latLngToWebMerc(double lat, double lng) {
  double x = 6378137.0 * lng * 0.017453292519943295;
  double a = lat * 0.017453292519943295;

  double y = 3189068.5 * log((1.0 + sin(a)) / (1.0 - sin(a)));
  return Point<T>(x, y);
}

// _____________________________________________________________________________
template <typename T>
inline Point<T> webMercToLatLng(double x, double y) {
  double lat = 114.591559026 * (atan(exp(y / 6378137.0)) - 0.78539825);
  double lon = x / 111319.4907932735677;
  return Point<T>(lon, lat);
}

// _____________________________________________________________________________
template <typename G1, typename G2>
inline double webMercMeterDist(const G1& a, const G2& b) {
  // euclidean distance on web mercator is in meters on equator,
  // and proportional to cos(lat) in both y directions

  double latA = 2 * atan(exp(a.getY() / 6378137.0)) - 1.5707965;
  double latB = 2 * atan(exp(b.getY() / 6378137.0)) - 1.5707965;

  return util::geo::dist(a, b) * cos((latA + latB) / 2.0);
}
}
}

#endif  // UTIL_GEO_GEO_H_
