// Copyright 2016, University of Freibur
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "log/Log.h"
#include "./PolyLine.h"
#include "../util/Geo.h"

using namespace transitmapper;
using namespace geo;
using util::geo::Point;
using util::geo::Line;
using util::geo::intersection;

// _____________________________________________________________________________
PolyLine::PolyLine() {

}

// _____________________________________________________________________________
PolyLine::PolyLine(const Point& from, const Point& to) {
  *this << from << to;
}

// _____________________________________________________________________________
PolyLine& PolyLine::operator<<(const Point& p) {
  _line.push_back(p);
  return *this;
}

// _____________________________________________________________________________
PolyLine& PolyLine::operator>>(const Point& p) {
  _line.insert(_line.begin(), p);
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


  if (_line.size() < 2) return;

  Line ret;
  Point lastP = _line.front();

  Point* lastIns = 0;
  Point* befLastIns = 0;

  for (size_t i = 1; i < _line.size(); i++) {
    Point curP = _line[i];

    double n1 = lastP.get<1>() - curP.get<1>();
    double n2 = curP.get<0>() - lastP.get<0>();
    double n = sqrt(n1*n1 + n2*n2);

    // if n == 0, the segment is effectively a point
    // we would get into all sorts of troubles if we tried to offset a point
    if (!(n > 0)) continue;

    n1 = n1 / n;
    n2 = n2 / n;

    lastP.set<0>(lastP.get<0>() + (n1 * units));
    lastP.set<1>(lastP.get<1>() + (n2 * units));

    curP.set<0>(curP.get<0>() + (n1 * units));
    curP.set<1>(curP.get<1>() + (n2 * units));

    if (lastIns && befLastIns &&
        util::geo::lineIntersects(*lastIns, *befLastIns, lastP, curP)) {
      *lastIns = util::geo::intersection(*lastIns, *befLastIns, lastP, curP);

      double d = util::geo::dist(lastP, *lastIns);
      double d2 = util::geo::distToSegment(*lastIns, *befLastIns, lastP);

      if (d > fabs(units) * 2 && d2 < d - (fabs(units))) {
        geo::PolyLine pl(*lastIns, *befLastIns);
        geo::PolyLine pll(*lastIns, curP);
        pl = pl.getSegment(0, (d - (fabs(units))) / pl.getLength());
        pll = pll.getSegment(0, (d - (fabs(units))) / pll.getLength());

        ret.push_back(pll.getLine().back());
        *lastIns = pl.getLine().back();

        ret.push_back(curP);
      } else {
        ret.push_back(curP);
      }
    } else {
      ret.push_back(lastP);
      ret.push_back(curP);
    }

    lastIns = &ret[ret.size() - 1];
    befLastIns = &ret[ret.size() - 2];

    lastP = _line[i];
  }

  _line = ret;

  simplify(1);
  fixTopology(fabs(2*3.14*units));
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

  // find a more performant way to clear the result of above
  ret.simplify(0);

  return ret;
}

// _____________________________________________________________________________
PointOnLine PolyLine::getPointAtDist(double atDist) const {
  if (atDist > getLength()) atDist = getLength();
  if (atDist < 0) atDist = 0;

  double dist = 0;

  if (_line.size() == 1) return PointOnLine(0, 0, _line[0]);

  const Point* last = &_line[0];

  for (size_t i = 1; i < _line.size(); i++) {
    const Point& cur = _line[i];
    double d = bgeo::distance(last, cur);
    dist += d;

    if (dist > atDist) {
      double p = (d - (dist - atDist));
      return PointOnLine(i-1, atDist / getLength(), interpolate(*last, cur, p));
    }

    last = &_line[i];
  }

  return PointOnLine(_line.size()-1, 1, _line.back());
}

// _____________________________________________________________________________
PointOnLine PolyLine::getPointAt(double at) const {
  at *= getLength();
  return getPointAtDist(at);
}

// _____________________________________________________________________________
Point PolyLine::interpolate(const Point& a,
  const Point& b, double p) const {
  double n1 = b.get<0>() - a.get<0>();
  double n2 = b.get<1>() - a.get<1>();
  double n = sqrt(n1*n1 + n2*n2);
  n1 = n1 / n;
  n2 = n2 / n;
  return Point(a.get<0>() + (n1 * p), a.get<1>() + (n2 * p));
}

// _____________________________________________________________________________
double PolyLine::distTo(const PolyLine& g) const {
  return bgeo::distance(_line, g.getLine());
}

// _____________________________________________________________________________
double PolyLine::distTo(const Point& p) const {
  return bgeo::distance(_line, p);
}

// _____________________________________________________________________________
double PolyLine::getLength() const {
  return bgeo::length(_line);
}

// _____________________________________________________________________________
PolyLine PolyLine::average(const std::vector<const PolyLine*>& lines,
    const std::vector<double>& weights) {
  bool weighted = lines.size() == weights.size();
  double stepSize;

  double longestLength = DBL_MIN;  // avoid recalc of length on each comparision
  for (const PolyLine* p : lines) {
    if (p->getLength() > longestLength) {
      longestLength = p->getLength();
    }
  }

  PolyLine ret;
  double total = 0;

  for (size_t i = 0; i < lines.size(); ++i) {
    if (weighted) {
      total += weights[i];
    } else {
      total += 1;
    }
  }

  stepSize = AVERAGING_STEP / longestLength;
  bool end = false;
  for (double a = 0; !end; a += stepSize) {
    if (a > 1) {
      a = 1;
      end = true;
    }
    double x = 0;
    double y = 0;

    for (size_t i = 0; i < lines.size(); ++i) {
      const PolyLine* pl = lines[i];
      Point p = pl->getPointAt(a).p;
      if (weighted) {
        x += p.get<0>() * weights[i];
        y += p.get<1>() * weights[i];
      } else {
        x += p.get<0>();
        y += p.get<1>();
      }
    }
    ret << Point(x / total, y / total);
  }

  ret.simplify(0);

  return ret;
}

// _____________________________________________________________________________
PolyLine PolyLine::average(const std::vector<const PolyLine*>& lines) {
  return average(lines, std::vector<double>());
}

// _____________________________________________________________________________
std::pair<size_t, double> PolyLine::nearestSegmentAfter(const Point& p,
    size_t a) const {
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
      totalDist += bgeo::distance(_line[i-2], _line[i-1]);
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

  if (totalLength > 0) {
    smallestDist /= totalLength;
  } else {
    smallestDist = 0;
  }

  return std::pair<size_t, double>(smallest, smallestDist);
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

  if (getLength() > 0) {
    bc.second += bgeo::distance(_line[bc.first], ret) / getLength();
  }

  return PointOnLine(bc.first, bc.second, ret);
}

// _____________________________________________________________________________
void PolyLine::simplify(double d) {
  util::geo::Line simpled;
  bgeo::simplify(_line, simpled, d);
  _line = simpled;
}

// _____________________________________________________________________________
void PolyLine::smoothenOutliers(double d) {
  if (_line.size() < 3) return;
  for (size_t i = 1; i < _line.size()-3; ++i) {
    double ang = util::geo::innerProd(_line[i], _line[i-1], _line[i+1]);

    if (util::geo::dist(_line[i], _line[i+1]) < d || util::geo::dist(_line[i], _line[i-1]) < d) {
      if (ang < 35) {
        _line.erase(_line.begin() + i);
      }
    }
  }
}

// _____________________________________________________________________________
bool PolyLine::equals(const PolyLine& rhs) const {
  // TODO: why 100? make global static or configurable or determine in some
  //       way!
  return equals(rhs, 100);
}

// _____________________________________________________________________________
bool PolyLine::operator==(const PolyLine& rhs) const {
  // TODO: why 100? make global static or configurable or determine in some
  //       way!
  return equals(rhs, 100);
}

// _____________________________________________________________________________
bool PolyLine::equals(const PolyLine& rhs, double dmax) const {
  // check if two lines are equal, THE DIRECTION DOES NOT MATTER HERE!!!!!

  if (_line.size() == 2 &&_line.size() == rhs.getLine().size()) {
    // trivial case, straight line, implement directly
    return (bgeo::distance(_line[0], rhs.getLine()[0]) < dmax &&
      bgeo::distance(_line.back(), rhs.getLine().back()) < dmax) ||
      (bgeo::distance(_line[0], rhs.getLine().back()) < dmax &&
      bgeo::distance(_line.back(), rhs.getLine()[0]) < dmax);
  } else {
    return contains(rhs, dmax) && rhs.contains(*this, dmax);
  }

  return true;
}

// _____________________________________________________________________________
bool PolyLine::contains(const PolyLine& rhs, double dmax) const {
  // check if two lines are equal, THE DIRECTION DOES NOT MATTER HERE!!!!!

  for (size_t i = 0; i < rhs.getLine().size(); ++i) {
    double d = bgeo::distance(rhs.getLine()[i], getLine());
    if (d > dmax) {
      return false;
    }
  }

  return true;
}

// _____________________________________________________________________________
void PolyLine::move(double vx, double vy) {
  for (size_t i = 0; i < _line.size(); i++) {
    _line[i].set<0>(_line[i].get<0>() + vx);
    _line[i].set<1>(_line[i].get<1>() + vy);
  }
}

// _____________________________________________________________________________
SharedSegments PolyLine::getSharedSegments(const PolyLine& pl, double dmax)
const {
  /**
   * Returns the segments this polyline share with pl
   * atm, this is a very simple distance-based algorithm
   *
   * TODO: use some mutation of frechet distance here..?
   */
  double STEP_SIZE = 2;
  double MAX_SKIPS = 4;
  double MIN_SEG_LENGTH = dmax / 2; // make this configurable!


  SharedSegments ret;

  if (distTo(pl) > dmax) return ret;

  bool in = false;
  double curDist = 0;
  double curTotalSegDist = 0;
  bool single = true;
  size_t skips;
  PointOnLine currentStartCand;
  PointOnLine currentEndCand;

  PointOnLine currentStartCandComp;
  PointOnLine currentEndCandComp;

  double comp = 0;
  double curSegDist = 0;
  for (size_t i = 1; i < _line.size(); ++i) {
    const Point& s = _line[i-1];
    const Point& e = _line[i];

    bool lastRound = false;

    double totalDist = bgeo::distance(s, e);
    while (curSegDist <= totalDist) {
      const Point& curPointer = interpolate(s, e, curSegDist);

      if (pl.distTo(curPointer) <= dmax) {
        PointOnLine curCompPointer = pl.projectOn(curPointer);
        PointOnLine curBackProjectedPointer = projectOn(curCompPointer.p);
        skips = 0;

        if (in) {
          currentEndCand = curBackProjectedPointer;
          currentEndCandComp = curCompPointer;

          single = false;

          comp = fabs(currentStartCand.totalPos * getLength() - currentEndCand.totalPos * getLength()) / fabs(currentStartCandComp.totalPos * pl.getLength() - currentEndCandComp.totalPos * pl.getLength());
        } else {
          in = true;
          currentStartCand = curBackProjectedPointer;
          currentStartCandComp = curCompPointer;
        }
      } else {
        if (in) {
          skips++;
          if (skips > MAX_SKIPS) { // TODO: make configurable
            if (comp > 0.8 && comp < 1.2 && !single &&
                (
                  fabs(currentStartCand.totalPos * getLength() - currentEndCand.totalPos * getLength()) > MIN_SEG_LENGTH &&
                  fabs(currentStartCandComp.totalPos * pl.getLength() - currentEndCandComp.totalPos * pl.getLength()) > MIN_SEG_LENGTH
                )
              ) {

              /**
              bool intersects =  false;
              if (util::geo::intersects(currentStartCand.p, currentEndCand.p, currentStartCandComp.p, currentEndCandComp.p)) {
                Point is = util::geo::intersection(currentStartCand.p, currentEndCand.p, currentStartCandComp.p, currentEndCandComp.p);
                double d1 = util::geo::dist(currentStartCand.p, is);
                double d2 = util::geo::dist(currentEndCand.p, is);

                if (d1 > dmax / 4 && d2 > dmax / 4) intersects = true;
              }
              **/

              ret.segments.push_back(SharedSegment(std::pair<PointOnLine, PointOnLine>(currentStartCand, currentStartCandComp),
                  std::pair<PointOnLine, PointOnLine>(currentEndCand, currentEndCandComp)));

              // TODO: only return the FIRST one, make this configuralbe
              return ret;
            }

            in = false;
            single = true;
          }
        }
      }

      if (curSegDist + STEP_SIZE > totalDist && !lastRound) {
        lastRound = true;
        double finalStep = totalDist - curSegDist - 0.0005;
        curSegDist += finalStep;
        curDist += finalStep;
      } else {
        curSegDist += STEP_SIZE;
        curDist += STEP_SIZE;
      }
    }

    curSegDist = curSegDist - totalDist;
    curTotalSegDist += totalDist;
  }

  if (comp > 0.8 && comp < 1.2 && in && !single &&
      (
        fabs(currentStartCand.totalPos * getLength() - currentEndCand.totalPos * getLength()) > MIN_SEG_LENGTH &&
        fabs(currentStartCandComp.totalPos * pl.getLength() - currentEndCandComp.totalPos * pl.getLength()) > MIN_SEG_LENGTH
      )
    ) {
    ret.segments.push_back(SharedSegment(std::pair<PointOnLine, PointOnLine>(currentStartCand, currentStartCandComp),
        std::pair<PointOnLine, PointOnLine>(currentEndCand, currentEndCandComp)));
  }

  return ret;
}

// _____________________________________________________________________________
std::set<PointOnLine, PointOnLineCompare> PolyLine::getIntersections(const PolyLine& g) const {
  std::set<PointOnLine, PointOnLineCompare> ret;

  for (size_t i = 1; i < g.getLine().size(); ++i) {
    // for each line segment, check if it intersects with a line segment in g
    const std::set<PointOnLine, PointOnLineCompare> a = getIntersections(g, i-1, i);
    ret.insert(a.begin(), a.end());
  }

  return ret;
}

// _____________________________________________________________________________
std::set<PointOnLine, PointOnLineCompare> PolyLine::getIntersections(const PolyLine& p,
    size_t a, size_t b) const {
  std::set<PointOnLine, PointOnLineCompare> ret;

  if (util::geo::dist(p.getLine()[a], p.getLine()[b]) == 0) {
    // we cannot intersect with a point
    return ret;
  }

  for (size_t i = 1; i < _line.size(); ++i) {
    if (util::geo::intersects(_line[i-1], _line[i], p.getLine()[a], p.getLine()[b])) {
      Point isect = intersection(_line[i-1], _line[i], p.getLine()[a], p.getLine()[b]);
      ret.insert(p.projectOn(isect));
    }
  }

  return ret;
}

// _____________________________________________________________________________
PolyLine PolyLine::getOrthoLineAtDist(double d, double length) const {
  Point avgP = getPointAtDist(d).p;

  double angle = util::geo::angBetween(getPointAtDist(d-5).p,
      getPointAtDist(d+5).p);

  double angleX1 = avgP.get<0>() + cos(angle + M_PI/2) * length/2;
  double angleY1 = avgP.get<1>() + sin(angle + M_PI/2) * length/2;

  double angleX2 = avgP.get<0>() + cos(angle + M_PI/2) * -length/2;
  double angleY2 = avgP.get<1>() + sin(angle + M_PI/2) * -length/2;

  geo::PolyLine pl(Point(angleX1, angleY1), Point(angleX2, angleY2));

  return pl;
}

// _____________________________________________________________________________
void PolyLine::empty() {
  _line.empty();
}

// _____________________________________________________________________________
std::pair<double, double> PolyLine::getSlopeBetween(double ad, double bd) const {
  PointOnLine a = getPointAt(ad);
  PointOnLine b = getPointAt(bd);

  double d = util::geo::dist(a.p, b.p);

  double dx = (b.p.get<0>() - a.p.get<0>()) / d;
  double dy = (b.p.get<1>() - a.p.get<1>()) / d;

  return std::pair<double, double>(dx, dy);
}

// _____________________________________________________________________________
std::pair<double, double> PolyLine::getSlopeBetweenDists(double ad, double bd)
const {
  return getSlopeBetween(ad / getLength(), bd / getLength());

}

// _____________________________________________________________________________
std::string PolyLine::getWKT() const {
  std::stringstream ss;
  ss << std::setprecision(12) << bgeo::wkt(_line);

  return ss.str();
}

// _____________________________________________________________________________
void PolyLine::fixTopology(double maxl) {
  double distA = 0;

  for (size_t i = 1; i < _line.size()-1; i++) {
    double distB = distA + util::geo::dist(_line[i-1], _line[i]) +
      util::geo::dist(_line[i], _line[i+1]);
    for (size_t j = i+2; j < _line.size(); j++) {
      if (util::geo::intersects(_line[i-1], _line[i], _line[j-1], _line[j])) {
        Point p = util::geo::intersection(_line[i-1], _line[i], _line[j-1], _line[j]);

        double posA = (util::geo::dist(_line[i-1], p) + distA);
        double posB = (util::geo::dist(_line[j-1], p) + distB);

        if (fabs(posA - posB) < maxl) {
          _line[i] = p;
          _line.erase(_line.begin() + i + 1, _line.begin() + j);
        }
      }

      distB += util::geo::dist(_line[j-1], _line[j]);
    }
    distA += util::geo::dist(_line[i-1], _line[i]);
  }
}

// _____________________________________________________________________________
void PolyLine::applyChaikinSmooth(size_t depth) {
  for (size_t i = 0; i < depth; i++) {
    Line smooth;

    smooth.push_back(_line.front());

    for (size_t i = 1; i < _line.size(); i++) {
      Point pA = _line[i-1];
      Point pB = _line[i];

      smooth.push_back(Point(0.75 * pA.get<0>() + 0.25 * pB.get<0>(), 0.75 * pA.get<1>() + 0.25 * pB.get<1>()));
      smooth.push_back(Point(0.25 * pA.get<0>() + 0.75 * pB.get<0>(), 0.25 * pA.get<1>() + 0.75 * pB.get<1>()));
    }

    smooth.push_back(_line.back());

    _line = smooth;
  }
}
