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

// TODO: maybe let this class inherit from a more generic geometry class
class PolyLine {

 public:
  PolyLine();
  PolyLine(const util::geo::Point& from, const util::geo::Point& to);

  PolyLine& operator<<(const util::geo::Point& p);

  void reverse();

  void offsetPerp(double units);

  PolyLine getPerpOffsetted(double units) const;

	const util::geo::Line& getLine() const;

  double distTo(const PolyLine& g) const;
  double distTo(const util::geo::Point& p) const;

  void getSharedSegments(const PolyLine& pl,
    std::vector<SharedSeg>* res) const;

  double getLength() const;

  // return point at dist
  util::geo::Point getPointAtDist(double dist) const;

  // return point at [0..1]
  util::geo::Point getPointAt(double dist) const;

  static PolyLine average(std::vector<const PolyLine*>& lines);

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
