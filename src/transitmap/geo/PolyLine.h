// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GEO_POLYLINE_H_
#define TRANSITMAP_GEO_POLYLINE_H_

#include <vector>
#include <string>
#include <ostream>
#include "../util/Geo.h"

namespace transitmapper {
namespace geo {

typedef std::pair<std::pair<size_t, size_t>, std::pair<size_t, size_t> > SharedSeg;

static const double MAX_EQ_DISTANCE = 15;
static const double AVERAGING_STEP = 20;

struct PointOnLine {
  PointOnLine(size_t i, double pos, const util::geo::Point& p)
  : lastIndex(i), totalPos(pos), p(p) {}
  size_t lastIndex;
  double totalPos;
  util::geo::Point p;
};

// TODO: maybe let this class inherit from a more generic geometry class
class PolyLine {

 public:
  PolyLine();
  PolyLine(const util::geo::Point& from, const util::geo::Point& to);

  PolyLine& operator<<(const util::geo::Point& p);

  void reverse();
  PolyLine getReversed() const;

  void offsetPerp(double units);
  void offsetPerp(double units, double xscale, double yscale);

  PolyLine getPerpOffsetted(double units) const;
  PolyLine getPerpOffsetted(double units, double xscale, double yscale) const;

	const util::geo::Line& getLine() const;

  double distTo(const PolyLine& g) const;
  double distTo(const util::geo::Point& p) const;

  void getSharedSegments(const PolyLine& pl,
    std::vector<SharedSeg>* res) const;

  double getLength() const;

  // return point at dist
  PointOnLine getPointAtDist(double dist) const;

  // return point at [0..1]
  PointOnLine getPointAt(double dist) const;

  PolyLine getSegment(double a, double b) const;
  PolyLine getSegment(const PointOnLine& start, const PointOnLine& end) const;
  PolyLine getSegment(const util::geo::Point& a, const util::geo::Point& b) const;

  static PolyLine average(std::vector<const PolyLine*>& lines);

  void simplify(double d);

	std::pair<size_t, double> nearestSegment(const util::geo::Point& p) const;
	std::pair<size_t, double> nearestSegmentAfter(const util::geo::Point& p, size_t after) const;

  PointOnLine projectOn(const util::geo::Point& p) const;
  PointOnLine projectOnAfter(const util::geo::Point& p, size_t after) const;

  // equality operator, will hold frechet-distance equality check in
  // the future
  bool operator==(const PolyLine& rhs) const;
 private:
  util::geo::Line _line;

  // TODO: maybe move this into util header
  util::geo::Point interpolate(const util::geo::Point& a, const util::geo::Point& b, double p) const;
};

}}

#endif  // TRANSITMAP_GEO_POLYLINE_H_
