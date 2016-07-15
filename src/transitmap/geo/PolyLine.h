// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GEO_POLYLINE_H_
#define TRANSITMAP_GEO_POLYLINE_H_

#include <string>
#include <ostream>
#include "../util/Geo.h"

namespace transitmapper {
namespace geo {

// TODO: maybe let this class inherit from a more generic geometry class
class PolyLine {

 public:
  PolyLine();

  PolyLine& operator<<(const util::geo::Point& p);

  void reverse();

  void offsetPerp(double units);

  PolyLine getPerpOffsetted(double units) const;

	const util::geo::Line& getLine() const;

 private:
  util::geo::Line _line;

};

}}

#endif  // TRANSITMAP_GEO_POLYLINE_H_
