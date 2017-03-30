// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GEO_POLYLINE_H_
#define TRANSITMAP_GEO_POLYLINE_H_

#include <ostream>
#include <string>
#include <vector>
#include "Geo.h"

namespace pbutil {
namespace geo {

static const double MAX_EQ_DISTANCE = 15;
static const double AVERAGING_STEP = 20;

struct LinePoint {
  LinePoint() : lastIndex(0), totalPos(-1), p() {}

  LinePoint(size_t i, double pos, const Point& p)
      : lastIndex(i), totalPos(pos), p(p) {}
  size_t lastIndex;
  double totalPos;
  Point p;
};

struct LinePointCmp {
  bool operator()(const LinePoint& lh, const LinePoint& rh) const {
    return lh.totalPos < rh.totalPos;
  }
};

typedef std::pair<LinePoint, LinePoint> LinePointPair;
typedef std::pair<LinePointPair, LinePointPair> SharedSegment;

struct SharedSegments {
  std::vector<SharedSegment> segments;
};

// TODO: maybe let this class inherit from a more generic geometry class
class PolyLine {
 public:
  PolyLine();
  PolyLine(const Point& from, const Point& to);

  PolyLine& operator<<(const Point& p);
  PolyLine& operator>>(const Point& p);

  void reverse();
  PolyLine getReversed() const;

  void offsetPerp(double units);

  PolyLine getPerpOffsetted(double units) const;

  const Line& getLine() const;

  double distTo(const PolyLine& g) const;
  double distTo(const Point& p) const;

  SharedSegments getSharedSegments(const PolyLine& pl, double dmax) const;

  double getLength() const;

  // return point at dist
  LinePoint getPointAtDist(double dist) const;

  // return point at [0..1]
  LinePoint getPointAt(double dist) const;

  PolyLine getSegment(double a, double b) const;
  PolyLine getSegmentAtDist(double dista, double distb) const;
  PolyLine getSegment(const LinePoint& start, const LinePoint& end) const;
  PolyLine getSegment(const Point& a, const Point& b) const;

  std::set<LinePoint, LinePointCmp> getIntersections(const PolyLine& g) const;

  static PolyLine average(const std::vector<const PolyLine*>& lines);
  static PolyLine average(const std::vector<const PolyLine*>& lines,
                          const std::vector<double>& weights);

  void simplify(double d);
  void empty();

  void smoothenOutliers(double d);

  std::pair<size_t, double> nearestSegment(const Point& p) const;
  std::pair<size_t, double> nearestSegmentAfter(const Point& p,
                                                size_t after) const;

  LinePoint projectOn(const Point& p) const;
  LinePoint projectOnAfter(const Point& p, size_t after) const;

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

  Point interpolate(const Point& a, const Point& b, double p) const;

  void fixTopology(double maxl);
  void applyChaikinSmooth(size_t depth);

 private:
  std::set<LinePoint, LinePointCmp> getIntersections(const PolyLine& p,
                                                     size_t a, size_t b) const;
  Line _line;
};
}  // namespace geo
}  // namespace pbutil

#endif  // TRANSITMAP_GEO_POLYLINE_H_
