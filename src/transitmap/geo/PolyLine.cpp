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
  // bound to go wrong at /some/ point (self intersections, numerical
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

// _____________________________________________________________________________
util::geo::Point PolyLine::getPointAtDist(double atDist) const {
  double dist = 0;

  if (_line.size() == 1) return _line[0];

  const util::geo::Point* last = &_line[0];

  for (size_t i = 1; i < _line.size(); i++) {
    const util::geo::Point& cur = _line[i];
    double d = boost::geometry::distance(last, cur);
    dist += d;

    if (dist > atDist) {
      double p = (d - (dist - atDist)) / d;
      return interpolate(*last, cur, p);
    }
    last = &_line[i];
  }

  return _line.back();
}

// _____________________________________________________________________________
util::geo::Point PolyLine::getPointAt(double atDist) const {
  atDist *= boost::geometry::length(_line);
  return getPointAtDist(atDist);
}

// _____________________________________________________________________________
util::geo::Point PolyLine::interpolate(const util::geo::Point& a,
  const util::geo::Point& b, double p) const {
  double n1 = b.get<0>() - a.get<0>();
  double n2 = b.get<1>() - a.get<1>();
  double n = sqrt(n1*n1 + n2*n2);
  n1 = n1 / n;
  n2 = n2 / n;
  return util::geo::Point(a.get<0>() + (n1 * p * n), a.get<1>() + (n2 * p * n));
}

// _____________________________________________________________________________
double PolyLine::distTo(const PolyLine& g) const {
  return boost::geometry::distance(_line, g.getLine());
}

// _____________________________________________________________________________
double PolyLine::distTo(const util::geo::Point& p) const {
  return boost::geometry::distance(_line, p);
}

// _____________________________________________________________________________
void PolyLine::getSharedSegments(const PolyLine& pl,
  std::vector<SharedSeg>* ret) const {

  double curTotalDist = 0;
  double curTotalDistNeg = 0;
  std::pair<size_t, size_t> curStartCand;
  std::pair<size_t, size_t> curEndCand;

  for (const util::geo::Point& p : _line) {
    // for now, only check anchor points, should be sufficient for most
    // lines

    // TODO


  }
}

// _____________________________________________________________________________
double PolyLine::getLength() const {
  return boost::geometry::length(_line);
}

// _____________________________________________________________________________
PolyLine PolyLine::average(std::vector<const PolyLine*>& lines) {
  double stepSize;
  const PolyLine* longest;
  double longestLength = 0;  // avoid recalc of length on each comparision
  for (const PolyLine* p : lines) {
    if (p->getLength() > longestLength) {
      longestLength = p->getLength();
      longest = p;
    }
  }

  PolyLine ret;
  stepSize = AVERAGING_STEP / longestLength;

  for (double a = 0; a <= 1.0; a += stepSize) {
    double x = 0;
    double y = 0;

    for (const PolyLine* pl : lines) {
      util::geo::Point p = pl->getPointAt(a);
      x += p.get<0>();
      y += p.get<1>();
    }
    ret << util::geo::Point(x / lines.size(), y / lines.size());
  }

  return ret;
}
