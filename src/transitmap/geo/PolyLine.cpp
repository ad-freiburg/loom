// Copyright 2016, University of Freibur
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
PolyLine::PolyLine(const util::geo::Point& from, const util::geo::Point& to) {
  *this << from << to;
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
PolyLine PolyLine::getReversed() const {
  PolyLine ret = *this;
  ret.reverse();
  return ret;
}

// _____________________________________________________________________________
const Line& PolyLine::getLine() const {
  return _line;
}

// _____________________________________________________________________________
PolyLine PolyLine::getPerpOffsetted(double units) const {
  return getPerpOffsetted(units, 1, 1);
}

// _____________________________________________________________________________
PolyLine PolyLine::getPerpOffsetted(double units, double xscale, double yscale) const {
  PolyLine p = *this;
  p.offsetPerp(units);
  return p;
}

// _____________________________________________________________________________
void PolyLine::offsetPerp(double units) {
  return offsetPerp(units, 1, 1);
}

// _____________________________________________________________________________
void PolyLine::offsetPerp(double units, double xscale, double yscale) {
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

    lastP.set<0>(lastP.get<0>() + (n1 * units / xscale));
    lastP.set<1>(lastP.get<1>() + (n2 * units / yscale));

    curP.set<0>(curP.get<0>() + (n1 * units / xscale));
    curP.set<1>(curP.get<1>() + (n2 * units / yscale));

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
PolyLine PolyLine::getSegment(double a, double b) const {
  if (a > b) a = b;
  PointOnLine start = getPointAt(a);
  PointOnLine end = getPointAt(b);

  return getSegment(start, end);
}

// _____________________________________________________________________________
PolyLine PolyLine::getSegment(const Point& a, const Point& b) const {
  PointOnLine start = projectOn(a);
  PointOnLine end = projectOnAfter(b, start.lastIndex);

  return getSegment(start, end);
}

// __________________________________________________________________
PolyLine PolyLine::getSegment(const PointOnLine& start, const PointOnLine& end)
const {
  PolyLine ret;

  ret << start.p;

  if (start.lastIndex+1 <= end.lastIndex) {
    ret._line.insert(
      ret._line.end(), _line.begin() + start.lastIndex+1, _line.begin() + end.lastIndex + 1
    );
  }
  ret << end.p;

  return ret;
}

// _____________________________________________________________________________
PointOnLine PolyLine::getPointAtDist(double atDist) const {
  double dist = 0;

  if (_line.size() == 1) return PointOnLine(0, 0, _line[0]);

  const util::geo::Point* last = &_line[0];

  for (size_t i = 1; i < _line.size(); i++) {
    const util::geo::Point& cur = _line[i];
    double d = boost::geometry::distance(last, cur);
    dist += d;

    if (dist > atDist) {
      double p = (d - (dist - atDist)) / d;
      return PointOnLine(i-1, atDist / getLength(), interpolate(*last, cur, p));
    }

    last = &_line[i];
  }

  return PointOnLine(_line.size()-1, 1, _line.back());
}

// _____________________________________________________________________________
PointOnLine PolyLine::getPointAt(double atDist) const {
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
double PolyLine::getLength() const {
  return boost::geometry::length(_line);
}

// _____________________________________________________________________________
PolyLine PolyLine::average(std::vector<const PolyLine*>& lines) {
  double stepSize;
  const PolyLine* longest = 0;
  double longestLength = DBL_MIN;  // avoid recalc of length on each comparision
  for (const PolyLine* p : lines) {
    if (p->getLength() > longestLength) {
      longestLength = p->getLength();
      longest = p;
    }
  }


  PolyLine ret;
  stepSize = AVERAGING_STEP / longestLength;
  bool end = false;
  for (double a = 0; !end; a += stepSize) {
    if (a > 1) {
      a = 1;
      end = true;
    }
    double x = 0;
    double y = 0;

    for (const PolyLine* pl : lines) {
      util::geo::Point p = pl->getPointAt(a).p;
      x += p.get<0>();
      y += p.get<1>();
    }
    ret << util::geo::Point(x / lines.size(), y / lines.size());
  }

  return ret;
}

// _____________________________________________________________________________
std::pair<size_t, double> PolyLine::nearestSegmentAfter(const Point& p, size_t a) const {
  // returns the index of the starting point of the nearest segment of p

  assert(a < _line.size());

  double totalLength = getLength();
  size_t smallest = a;
  double totalDist = 0;
  double dist = DBL_MAX;
  double smallestDist = 0;

  for (size_t i = smallest + 1; i < _line.size(); i++) {
    Point startP(_line[i-1]);
    Point endP(_line[i]);

    if (i > 1) {
      totalDist += boost::geometry::distance(_line[i-2], _line[i-1]);
    }

    double curDist = util::geo::distToSegment(
      startP,
      endP,
      p
    );

    if (curDist < dist) {
      dist = curDist;
      smallest = i-1;
      smallestDist = totalDist;
    }
  }

  return std::pair<size_t, double>(smallest, smallestDist / totalLength);
}

// _____________________________________________________________________________
std::pair<size_t, double> PolyLine::nearestSegment(const Point& p) const {
  return nearestSegmentAfter(p, 0);
}

// _____________________________________________________________________________
PointOnLine PolyLine::projectOn(const Point& p) const {
  return projectOnAfter(p, 0);
}

// _____________________________________________________________________________
PointOnLine PolyLine::projectOnAfter(const Point& p, size_t a) const {
  assert(a < _line.size());
  std::pair<size_t, double> bc = nearestSegmentAfter(p, a);

  Point ret = util::geo::projectOn(
    _line[bc.first],
    p,
    _line[bc.first+1]
  );

  bc.second += boost::geometry::distance(_line[bc.first], ret) / getLength();

  return PointOnLine(bc.first, bc.second, ret);
}

// _____________________________________________________________________________
void PolyLine::simplify(double d) {
  util::geo::Line simpled;
  boost::geometry::simplify(_line, simpled, d);
  _line = simpled;
}

// _____________________________________________________________________________
bool PolyLine::operator==(const PolyLine& rhs) const {
  // check if two lines are equal, THE DIRECTION DOES NOT MATTER HERE!!!!!

  // TODO: why 100? make global static or configurable or determine in some
  //       way!
  double DMAX = 100;

  if (_line.size() == 2 &&_line.size() == rhs.getLine().size()) {
    // trivial case, straight line, implement directly
    return (boost::geometry::distance(_line[0], rhs.getLine()[0]) < DMAX &&
      boost::geometry::distance(_line.back(), rhs.getLine().back()) < DMAX) ||
      (boost::geometry::distance(_line[0], rhs.getLine().back()) < DMAX &&
      boost::geometry::distance(_line.back(), rhs.getLine()[0]) < DMAX);
  } else {
    return contains(rhs) && rhs.contains(*this);
  }

  return true;
}

// _____________________________________________________________________________
bool PolyLine::contains(const PolyLine& rhs) const {
  // check if two lines are equal, THE DIRECTION DOES NOT MATTER HERE!!!!!

  // TODO: why 100? make global static or configurable or determine in some
  //       way!
  double DMAX = 100;

  for (size_t i = 0; i < rhs.getLine().size(); ++i) {
    double d = boost::geometry::distance(rhs.getLine()[i], getLine());
    if (d > DMAX) {
      return false;
    }
  }

  return true;
}

// _____________________________________________________________________________ 
SharedSegments PolyLine::getSharedSegments(const PolyLine& pl) const {
  /**
   * Returns the segments this polyline share with pl
   * atm, this is a very simple distance-based algorithm
   *
   * TODO: use some mutation of frechet distance here
   */
  double STEP_SIZE = 20;
  double DMAX = 50;
  SharedSegments ret;

  bool in = false;
  bool notSingle = false;
  double segDist = 0;
  double curDist = 0;
  double curTotalSegDist = 0;
  size_t skips;
  PointOnLine currentStartCand;
  PointOnLine currentEndCand;

  double curSegDist = 0;
  for (size_t i = 1; i < _line.size(); ++i) {
    const Point& s = _line[i-1];
    const Point& e = _line[i];

    double totalDist = boost::geometry::distance(s, e);
    while (curSegDist < totalDist) {
      const Point& curPointer = interpolate(s, e, curSegDist / totalDist);
      auto nearestSeg = pl.nearestSegment(curPointer);


      // curPointer is the current point we are checking on THIS polyline
      // nearestSeg is the nearest total segment on pl we are checking against
      // we now have to check whether curPointer "contained"  in nearestSegment

      double ang = util::geo::innerProd(pl.getLine()[nearestSeg.first], curPointer, pl.getLine()[nearestSeg.first + 1]);
      double ang2 = util::geo::innerProd(pl.getLine()[nearestSeg.first + 1], curPointer, pl.getLine()[nearestSeg.first]);

      if ((ang) <= 90 && (ang2) <= 90 &&
          pl.distTo(curPointer) <= DMAX) {
        skips = 0;
        if (in) {
          // update currendEndCand
          if (notSingle) segDist += util::geo::dist(currentEndCand.p, curPointer);
          else segDist += util::geo::dist(currentStartCand.p, curPointer);
          notSingle = true;
          currentEndCand = PointOnLine(i-1, (curTotalSegDist + curSegDist) / getLength(), curPointer);
        } else {
          in = true;
          currentStartCand = PointOnLine(i-1, (curTotalSegDist + curSegDist) / getLength(), curPointer);
        }
      } else {
        if (in) {
          skips++;
          if (skips > 5) {
            if (notSingle && segDist > DMAX) {
              ret.segments.push_back(std::pair<PointOnLine, PointOnLine>(currentStartCand, currentEndCand));

              // TODO: only return the FIRST one, make this configuralbe
              return ret;
            } else {
            }
            in = false;
            notSingle = false;
            segDist = 0;
          }
        } else {
          // do nothing
        }
      }

      curSegDist += STEP_SIZE;
      curDist += STEP_SIZE;

    }

    curSegDist = curSegDist - totalDist;
    curTotalSegDist += totalDist;
  }

  if (in && notSingle && segDist > DMAX) {
    ret.segments.push_back(std::pair<PointOnLine, PointOnLine>(currentStartCand, currentEndCand));
  }
  return ret;
}

