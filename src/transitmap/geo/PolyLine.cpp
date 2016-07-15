// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include "PolyLine.h"

using namespace transitmapper;
using namespace geo;
using namespace util::geo;

// _____________________________________________________________________________
PolyLine::PolyLine() {

}

// _____________________________________________________________________________
PolyLine& PolyLine::operator<<(const util::geo::Point& p) {
  _line.push_back(p);
  return *this;
}

// _____________________________________________________________________________
void PolyLine::reverse() {
  std::reverse(_line.begin(), _line.end());
}

// _____________________________________________________________________________
const Line& PolyLine::getLine() const {
	return _line;
}

// _____________________________________________________________________________
PolyLine PolyLine::getPerpOffsetted(double units) const {
  PolyLine p = *this;
  p.offsetPerp(units);
  return p;
}

// _____________________________________________________________________________
void PolyLine::offsetPerp(double units) {
  // calculate perpendicular offset of a polyline
  //
  // there doesn't seem to be any library which reliably does that,
  // so we do it ourself here until we find one...
  // boost::geometry only supports buffering a line, resulting in a
  // polygon. An offsetted line is part of that polygon, but retrieving
  // it reliably could result in some geometrical diffing hocus pocus which is
  // bound to go wrong at _some_ point (self intersections, numerical
  // instability etc)
	Line ret;

	if (_line.size() < 2) return; // TODO: throw error

	Point lastP = _line.front();

  Point* lastIns = 0;
  Point* befLastIns = 0;

	for (size_t i = 1; i < _line.size(); i++) {
	  Point curP = _line[i];

		double n1 = lastP.get<1>() - curP.get<1>();
		double n2 = curP.get<0>() - lastP.get<0>();
		double n = sqrt(n1*n1 + n2*n2);
		n1 = n1 / n;
		n2 = n2 / n;

    lastP.set<0>(lastP.get<0>() + (n1 * units));
		lastP.set<1>(lastP.get<1>() + (n2 * units));

		curP.set<0>(curP.get<0>() + (n1 * units));
		curP.set<1>(curP.get<1>() + (n2 * units));

    if (lastIns && befLastIns &&
        lineIntersects(*lastIns, *befLastIns, lastP, curP)) {
      *lastIns = intersection(*lastIns, *befLastIns, lastP, curP);
      ret.push_back(curP);
    } else {
      ret.push_back(lastP);
      ret.push_back(curP);
    }

    lastIns = &ret[ret.size() - 1];
    befLastIns = &ret[ret.size() - 2];

		lastP = _line[i];
	}

	_line = ret;
}
