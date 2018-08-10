// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_POLYLINE_H_
#define UTIL_GEO_POLYLINE_H_

#include <cfloat>
#include <ostream>
#include <iomanip>
#include <string>
#include <set>
#include <vector>
#include "Geo.h"

namespace util {
namespace geo {

static const double MAX_EQ_DISTANCE = 15;
static const double AVERAGING_STEP = 20;

// legacy code, will be removed in the future

template <typename T>
struct LinePoint {
  LinePoint() : lastIndex(0), totalPos(-1), p() {}

  LinePoint(size_t i, double pos, const Point<T>& p)
      : lastIndex(i), totalPos(pos), p(p) {}
  size_t lastIndex;
  double totalPos;
  Point<T> p;
};

template <typename T>
struct LinePointCmp {
  bool operator()(const LinePoint<T>& lh, const LinePoint<T>& rh) const {
    return lh.totalPos < rh.totalPos;
  }
};


template <typename T>
using LinePointPair = std::pair<LinePoint<T>, LinePoint<T>>;
template <typename T>
using SharedSegment = std::pair<LinePointPair<T>, LinePointPair<T>>;

template <typename T>
struct SharedSegments {
  std::vector<SharedSegment<T>> segments;
};

template <typename T>
class PolyLine {
 public:
  PolyLine();
  PolyLine(const Point<T>& from, const Point<T>& to);
  PolyLine(const Line<T>& l);

  PolyLine& operator<<(const Point<T>& p);
  PolyLine& operator>>(const Point<T>& p);

  void reverse();
  PolyLine getReversed() const;

  void offsetPerp(double units);

  PolyLine getPerpOffsetted(double units) const;

  const Line<T>& getLine() const;

  double distTo(const PolyLine<T>& g) const;
  double distTo(const Point<T>& p) const;

  SharedSegments<T> getSharedSegments(const PolyLine<T>& pl, double dmax) const;

  double getLength() const;

  // return point at dist
  LinePoint<T> getPointAtDist(double dist) const;

  // return point at [0..1]
  LinePoint<T> getPointAt(double dist) const;

  PolyLine<T> getSegment(double a, double b) const;
  PolyLine<T> getSegmentAtDist(double dista, double distb) const;
  PolyLine<T> getSegment(const LinePoint<T>& start, const LinePoint<T>& end) const;
  PolyLine<T> getSegment(const Point<T>& a, const Point<T>& b) const;

  std::set<LinePoint<T>, LinePointCmp<T>> getIntersections(const PolyLine<T>& g) const;

  static PolyLine<T> average(const std::vector<const PolyLine<T>*>& lines);
  static PolyLine<T> average(const std::vector<const PolyLine<T>*>& lines,
                          const std::vector<double>& weights);

  void simplify(double d);
  void empty();

  void smoothenOutliers(double d);

  std::pair<size_t, double> nearestSegment(const Point<T>& p) const;
  std::pair<size_t, double> nearestSegmentAfter(const Point<T>& p,
                                                size_t after) const;

  LinePoint<T> projectOn(const Point<T>& p) const;
  LinePoint<T> projectOnAfter(const Point<T>& p, size_t after) const;

  void move(double vx, double vy);

  std::pair<double, double> getSlopeBetween(double ad, double bd) const;
  std::pair<double, double> getSlopeBetweenDists(double ad, double bd) const;

  // equality operator, will hold frechet-distance equality check in
  // the dmax
  bool operator==(const PolyLine& rhs) const;
  bool contains(const PolyLine& rhs, double dmax) const;
  bool equals(const PolyLine& rhs) const;
  bool equals(const PolyLine& rhs, double dmax) const;

  std::string getWKT() const;

  PolyLine getOrthoLineAtDist(double d, double lengt) const;

  Point<T> interpolate(const Point<T>& a, const Point<T>& b, double p) const;

  void fixTopology(double maxl);
  void applyChaikinSmooth(size_t depth);

  const Point<T>& front() const;
  const Point<T>& back() const;

 private:
  std::set<LinePoint<T>, LinePointCmp<T>> getIntersections(const PolyLine& p,
                                                     size_t a, size_t b) const;
  Line<T> _line;
};

#include "util/geo/PolyLine.tpp"

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_POLYLINE_H_
