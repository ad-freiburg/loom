// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include "../graph/OrderingConfiguration.h"
#include "./Edge.h"
#include "./Node.h"
#include "./Route.h"
#include "./TransitGraph.h"
#include "pbutil/geo/BezierCurve.h"
#include "pbutil/geo/Geo.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;

using pbutil::geo::Point;
using pbutil::geo::BezierCurve;
using pbutil::geo::PolyLine;
using pbutil::geo::Line;

// _____________________________________________________________________________
Point NodeFront::getTripOccPos(const Route* r,
                               const graph::Configuration& c) const {
  return getTripOccPos(r, c, false);
}

// _____________________________________________________________________________
Point NodeFront::getTripOccPos(const Route* r,
                               const graph::Configuration& c,
                               bool origGeom) const {
  RouteOccWithPos rop;

  assert(c.find(edge) != c.end());

  rop = edge->getRouteOccWithPosUnder(r, c.find(edge)->second);

  if (rop.first) {
    return getTripPos(edge, rop.second, n == edge->getTo(), origGeom);
  }

  throw std::runtime_error("route does not occur in this edge");
}

// _____________________________________________________________________________
Point NodeFront::getTripPos(const Edge* e, size_t pos, bool inv) const {
  return getTripPos(e, pos, inv, false);
}

// _____________________________________________________________________________
Point NodeFront::getTripPos(const Edge* e, size_t pos, bool inv, bool origG) const {
  double p;
  if (!inv) {
    p = (e->getWidth() + e->getSpacing()) * pos + e->getWidth() / 2;
  } else {
    p = (e->getWidth() + e->getSpacing()) * (e->getCardinality() - 1 - pos) +
        e->getWidth() / 2;
  }

  // use interpolate here directly for speed
  if (origG) {
    return origGeom.interpolate(origGeom.getLine().front(), origGeom.getLine().back(), p);
  } else {
    return geom.interpolate(geom.getLine().front(), geom.getLine().back(), p);
  }
}

// _____________________________________________________________________________
double Node::getMaxNodeFrontWidth() const {
  double ret = 0;
  for (const NodeFront& g : _mainDirs) {
    if (g.edge->getTotalWidth() > ret) ret = g.edge->getTotalWidth();
  }
  return ret;
}

// _____________________________________________________________________________
size_t Node::getMaxNodeFrontCardinality() const {
  size_t ret = 0;
  for (const NodeFront& g : _mainDirs) {
    if (g.edge->getCardinality() > ret) ret = g.edge->getCardinality();
  }
  return ret;
}

// _____________________________________________________________________________
Node::Node(const std::string& id, Point pos) : _id(id), _pos(pos) {}

// _____________________________________________________________________________
Node::Node(const std::string& id, double x, double y) : _id(id), _pos(x, y) {}

// _____________________________________________________________________________
Node::Node(const std::string& id, Point pos, StationInfo s)
    : _id(id), _pos(pos) {
  addStop(s);
}

// _____________________________________________________________________________
Node::Node(const std::string& id, double x, double y, StationInfo s)
    : _id(id), _pos(x, y) {
  addStop(s);
}

// _____________________________________________________________________________
Node::~Node() {
  for (auto e = _adjListOut.begin(); e != _adjListOut.end();) {
    Edge* eP = *e;

    if (eP->getFrom() == this) {
      // careful with invalidating iterators
      e = _adjListOut.erase(e);
    } else {
      eP->getFrom()->removeEdge(eP);
      e++;
    }

    eP->getTo()->removeEdge(eP);

    delete eP;
  }

  for (auto e = _adjListIn.begin(); e != _adjListIn.end();) {
    Edge* eP = *e;

    if (eP->getTo() == this) {
      // careful with invalidating iterators
      e = _adjListIn.erase(e);
    } else {
      eP->getTo()->removeEdge(eP);
      e++;
    }

    eP->getFrom()->removeEdge(eP);

    delete eP;
  }
}

// _____________________________________________________________________________
void Node::addStop(StationInfo s) { _stops.push_back(s); }

// _____________________________________________________________________________
const std::vector<StationInfo>& Node::getStops() const { return _stops; }

// _____________________________________________________________________________
void Node::addEdge(Edge* e) {
  if (e->getFrom() == this) _adjListOut.insert(e);
  if (e->getTo() == this) _adjListIn.insert(e);
}

// _____________________________________________________________________________
void Node::removeEdge(Edge* e) {
  if (e->getFrom() == this) _adjListOut.erase(e);
  if (e->getTo() == this) _adjListIn.erase(e);

  for (size_t i = 0; i < _mainDirs.size(); i++) {
    if (_mainDirs[i].edge == e) {
      _mainDirs.erase(_mainDirs.begin() + i);
    }
  }

  // TODO: remove from _routeConnExceptions
}

// _____________________________________________________________________________
const Point& Node::getPos() const { return _pos; }

// _____________________________________________________________________________
void Node::setPos(const Point& p) { _pos = p; }

// _____________________________________________________________________________
const std::string& Node::getId() const { return _id; }

// _____________________________________________________________________________
void Node::addMainDir(NodeFront f) { _mainDirs.push_back(f); }

// _____________________________________________________________________________
const NodeFront* Node::getNodeFrontFor(const Edge* e) const {
  for (auto& nf : getMainDirs()) {
    if (nf.edge == e) {
      return &nf;
    }
  }

  return 0;
}

// _____________________________________________________________________________
double Node::getScore(double inStatPen, double sameSegCrossPen,
                      double diffSegCrossPen, double splittingPen,
                      bool crossAdjPen, bool splitAdjPen,
                      const graph::Configuration& cfg) const {
  return getCrossingScore(cfg, inStatPen, sameSegCrossPen, diffSegCrossPen,
                          crossAdjPen) +
         getSeparationScore(cfg, inStatPen, splittingPen, splitAdjPen);
}

// _____________________________________________________________________________
size_t Node::getNumCrossings(const Configuration& c) const {
  return getCrossingScore(c, 1, 1, 1, false);
}

// _____________________________________________________________________________
double Node::getCrossingScore(const Configuration& c, double inStatPen,
                              double sameSegPen, double diffSegPen,
                              bool adjpen) const {
  std::vector<InnerGeometry> igs = getInnerGeometries(c, -1);
  size_t ret = 0;

  for (size_t i = 0; i < igs.size(); ++i) {
    for (size_t j = i + 1; j < igs.size(); ++j) {
      const InnerGeometry& iga = igs[i];
      const InnerGeometry& igb = igs[j];

      if (iga.from.front == 0 || iga.to.front == 0 || igb.from.front == 0 || igb.to.front == 0) continue;

      if (iga.from.front == igb.from.front && iga.slotFrom == igb.slotFrom)
        continue;
      if (iga.from.front == igb.to.front && iga.slotFrom == igb.slotTo)
        continue;
      if (iga.to.front == igb.to.front && iga.slotTo == igb.slotTo) continue;
      if (iga.to.front == igb.from.front && iga.slotTo == igb.slotFrom)
        continue;

      bool sameSeg =
          (iga.from.front == igb.from.front && iga.to.front == igb.to.front) ||
          (iga.to.front == igb.from.front && iga.from.front == igb.to.front);

      if (pbutil::geo::intersects(
              iga.geom.getLine().front(), iga.geom.getLine().back(),
              igb.geom.getLine().front(), igb.geom.getLine().back()) ||
          bgeo::distance(iga.geom.getLine(), igb.geom.getLine()) < 1) {
        ret += 1 * getCrossingPenalty(
                       inStatPen, sameSeg ? sameSegPen : diffSegPen, adjpen);
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
size_t Node::getNumSeparations(const Configuration& c) const {
  return getSeparationScore(c, 1, 1, false);
}

// _____________________________________________________________________________
size_t Node::getSeparationScore(const Configuration& c, double inStatPen,
                                double pen, bool adjpen) const {
  size_t ret = 0;
  for (auto nf : getMainDirs()) {
    const Edge* e = nf.edge;
    std::vector<std::pair<const Route*, const Route*> > curPairs;
    for (size_t i = 0; i < c.find(e)->second.size() - 1; i++) {
      size_t p = c.find(e)->second[i];
      curPairs.push_back(std::pair<const Route*, const Route*>(
          e->getTripsUnordered().at(p).route,
          e->getTripsUnordered().at(c.find(e)->second[i + 1]).route));
    }

    for (auto p : curPairs) {
      for (auto nf : getMainDirs()) {
        const Edge* f = nf.edge;
        if (e == f) continue;

        if (f->containsRoute(p.first) && f->containsRoute(p.second) &&
            connOccurs(p.first, e, f) && connOccurs(p.second, e, f)) {
          if (abs(int(f->getRouteOccWithPosUnder(p.first, c.find(f)->second)
                          .second) -
                  int(f->getRouteOccWithPosUnder(p.second, c.find(f)->second)
                          .second)) > 1) {
            ret += 1 * getSplittingPenalty(inStatPen, pen, adjpen);
          }
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::set<Edge*> Node::getAdjList() const {
  std::set<Edge*> ret;
  ret.insert(getAdjListIn().begin(), getAdjListIn().end());
  ret.insert(getAdjListOut().begin(), getAdjListOut().end());

  return ret;
}

// _____________________________________________________________________________
std::vector<Partner> Node::getPartners(const NodeFront* f,
                                       const RouteOccurance& ro) const {
  std::vector<Partner> ret;
  for (const auto& nf : getMainDirs()) {
    if (&nf == f) continue;

    for (const RouteOccurance& to :
         nf.edge->getContinuedRoutesIn(this, ro.route, ro.direction, f->edge)) {
      Partner p(f, nf.edge, to.route);
      p.front = &nf;
      p.edge = nf.edge;
      p.route = to.route;
      ret.push_back(p);
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<InnerGeometry> Node::getInnerGeometries(const Configuration& c,
                                                    double prec) const {
  std::vector<InnerGeometry> ret;
  std::map<const Route*, std::set<const NodeFront*> > processed;

  for (size_t i = 0; i < getMainDirs().size(); ++i) {
    const graph::NodeFront& nf = getMainDirs()[i];

    const std::vector<size_t>* ordering = &c.find(nf.edge)->second;

    for (size_t j : *ordering) {
      const RouteOccurance& routeOcc = (*nf.edge->getTripsUnordered())[j];
      Partner o(&nf, nf.edge, routeOcc.route);

      std::vector<graph::Partner> partners = getPartners(&nf, routeOcc);

      for (const graph::Partner& p : partners) {
        if (processed[routeOcc.route].find(p.front) !=
            processed[routeOcc.route].end()) {
          continue;
        }

        auto is = getInnerStraightLine(c, o, p);
        if (is.geom.getLength() == 0) continue;

        if (prec > 0) {
          ret.push_back(getInnerBezier(c, o, p, prec));
        } else {
          ret.push_back(is);
        }
      }

      // handle lines where this node is the terminus
      if (partners.size() == 0) {
        auto is = getTerminusStraightLine(c, o);
        if (is.geom.getLength() > 0) {
          if (prec > 0) {
            ret.push_back(getTerminusBezier(c, o, prec));
          } else {
            ret.push_back(is);
          }
        }
      }

      processed[routeOcc.route].insert(&nf);
    }
  }

  return ret;
}

// _____________________________________________________________________________
InnerGeometry Node::getTerminusStraightLine(
    const Configuration& c, const graph::Partner& partnerFrom) const {
  Point p = partnerFrom.front->getTripOccPos(partnerFrom.route, c, false);
  Point pp = partnerFrom.front->getTripOccPos(partnerFrom.route, c, true);

  size_t s = partnerFrom.edge
                 ->getRouteOccWithPosUnder(partnerFrom.route,
                                           c.find(partnerFrom.edge)->second)
                 .second;
  size_t ss = partnerFrom.edge
                  ->getRouteOccWithPosUnder(partnerFrom.route,
                                            c.find(partnerFrom.edge)->second)
                  .second;

  return InnerGeometry(PolyLine(p, pp), partnerFrom, Partner(), s, ss);
}

// _____________________________________________________________________________
InnerGeometry Node::getInnerStraightLine(
    const Configuration& c, const graph::Partner& partnerFrom,
    const graph::Partner& partnerTo) const {
  Point p = partnerFrom.front->getTripOccPos(partnerFrom.route, c);
  Point pp = partnerTo.front->getTripOccPos(partnerTo.route, c);

  size_t s = partnerFrom.edge
                 ->getRouteOccWithPosUnder(partnerFrom.route,
                                           c.find(partnerFrom.edge)->second)
                 .second;
  size_t ss = partnerTo.edge
                  ->getRouteOccWithPosUnder(partnerTo.route,
                                            c.find(partnerTo.edge)->second)
                  .second;

  return InnerGeometry(PolyLine(p, pp), partnerFrom, partnerTo, s, ss);
}

// _____________________________________________________________________________
InnerGeometry Node::getTerminusBezier(const Configuration& cf,
                                   const graph::Partner& partnerFrom,
                                   double prec) const {
  InnerGeometry ret = getTerminusStraightLine(cf, partnerFrom);
  Point p = ret.geom.getLine().front();
  Point pp = ret.geom.getLine().back();
  double d = pbutil::geo::dist(p, pp) / 2;

  Point b = p;
  Point c = pp;
  std::pair<double, double> slopeA, slopeB;

  assert(partnerFrom.front->edge->getGeom().getLength() > 5);

  if (partnerFrom.front->edge->getTo() == this) {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(
        partnerFrom.front->edge->getGeom().getLength() - 5,
        partnerFrom.front->edge->getGeom().getLength());
  } else {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(5, 0);
  }

  b = Point(p.get<0>() + slopeA.first * d, p.get<1>() + slopeA.second * d);
  c = Point(pp.get<0>(), pp.get<1>());

  BezierCurve bc(p, b, c, pp);
  ret.geom = bc.render(prec);

  return ret;
}

// _____________________________________________________________________________
InnerGeometry Node::getInnerBezier(const Configuration& cf,
                                   const graph::Partner& partnerFrom,
                                   const graph::Partner& partnerTo,
                                   double prec) const {
  InnerGeometry ret = getInnerStraightLine(cf, partnerFrom, partnerTo);
  Point p = ret.geom.getLine().front();
  Point pp = ret.geom.getLine().back();
  double d = pbutil::geo::dist(p, pp);

  Point b = p;
  Point c = pp;
  std::pair<double, double> slopeA, slopeB;

  assert(partnerFrom.front->edge->getGeom().getLength() > 5);
  assert(partnerTo.front->edge->getGeom().getLength() > 5);

  if (partnerFrom.front->edge->getTo() == this) {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(
        partnerFrom.front->edge->getGeom().getLength() - 5,
        partnerFrom.front->edge->getGeom().getLength());
  } else {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(5, 0);
  }

  if (partnerTo.front->edge->getTo() == this) {
    slopeB = partnerTo.front->edge->getGeom().getSlopeBetweenDists(
        partnerTo.front->edge->getGeom().getLength() - 5,
        partnerTo.front->edge->getGeom().getLength());
  } else {
    slopeB = partnerTo.front->edge->getGeom().getSlopeBetweenDists(5, 0);
  }

  double da = 1;
  double db = 1;

  Point pa = Point(p.get<0>() - slopeA.first * d * 2, p.get<1>() - slopeA.second * d * 2);
  Point pb = Point(pp.get<0>() - slopeB.first * d * 2, pp.get<1>() - slopeB.second * d * 2);
  Point ppa = Point(p.get<0>() + slopeA.first * d * 2, p.get<1>() + slopeA.second * d * 2);
  Point ppb = Point(pp.get<0>() + slopeB.first * d * 2, pp.get<1>() + slopeB.second * d * 2);


  if (d > 0.001 && pbutil::geo::intersects(pa, ppa, pb, ppb)) {
    Point isect = pbutil::geo::intersection(pa, ppa, pb, ppb);

    if (!std::isnan(isect.get<0>()) && !std::isnan(isect.get<1>())) {
      double degAng = pbutil::geo::innerProd(isect, pa, ppb);
      double ang = cos(degAng / (180 / M_PI));

      ang = ang * ang;

      if (std::isnan(ang)) ang = 1;

      double dar = pbutil::geo::dist(isect, p);
      double dbr = pbutil::geo::dist(isect, pp);

      double avg = (dar + dbr) / 2;

      da = (dar * 2 * (1 - ang) + avg * 2 * (ang)) / 2;
      db = (dbr * 2 * (1 - ang) + avg * 2 * (ang)) / 2;
    }
  }

  if ((da + db) < 0.001 || (da + db) > d * 2) {
    da = 1;
    db = 1;
  }

  b = Point(p.get<0>() + slopeA.first * d * (da / (da+db)), p.get<1>() + slopeA.second * d * (da / (da+db)));
  c = Point(pp.get<0>() + slopeB.first * d * (db / (da+db)), pp.get<1>() + slopeB.second * d * (db / (da+db)));

  BezierCurve bc(p, b, c, pp);
  ret.geom = bc.render(prec);

  return ret;
}

// _____________________________________________________________________________
size_t Node::getConnCardinality() const {
  size_t ret = 0;
  std::map<const Route*, std::set<const NodeFront*> > processed;

  for (size_t i = 0; i < getMainDirs().size(); ++i) {
    const graph::NodeFront& nf = getMainDirs()[i];

    for (size_t j = 0; j < nf.edge->getCardinality(); j++) {
      const RouteOccurance& routeOcc = (*nf.edge->getTripsUnordered())[j];

      std::vector<graph::Partner> partners = getPartners(&nf, routeOcc);

      for (const graph::Partner& p : partners) {
        if (processed[routeOcc.route].find(p.front) !=
            processed[routeOcc.route].end()) {
          continue;
        }
        ret++;
      }

      processed[routeOcc.route].insert(&nf);
    }
  }

  return ret;
}

// _____________________________________________________________________________
void Node::generateStationHull(double d, bool useSimple) {
  if (getMainDirs().size() == 0) return;
  _stationHull = getConvexFrontHull(d, true, useSimple);
}

// _____________________________________________________________________________
Polygon Node::getStationHull() const { return _stationHull; }

// _____________________________________________________________________________
Polygon Node::getConvexFrontHull(double d, bool rectangulize,
    bool simpleRenderForTwoEdgeNodes) const {
  double cd = d;

  MultiPolygon ret;
  double pointsPerCircle = 36;
  bgeo::strategy::buffer::distance_symmetric<double> distanceStrat(d);
  bgeo::strategy::buffer::join_round joinStrat(pointsPerCircle);
  bgeo::strategy::buffer::end_round endStrat(pointsPerCircle);
  bgeo::strategy::buffer::point_circle circleStrat(pointsPerCircle);
  bgeo::strategy::buffer::side_straight sideStrat;

  if (!simpleRenderForTwoEdgeNodes || getMainDirs().size() != 2) {
    MultiLine l;
    for (auto& nf : getMainDirs()) {
      l.push_back(
          nf.geom
              .getSegment((cd / 2) / nf.geom.getLength(),
                          (nf.geom.getLength() - cd / 2) / nf.geom.getLength())
              .getLine());
    }

    Polygon hull;
    bgeo::convex_hull(l, hull);

    if (rectangulize && getMaxNodeFrontCardinality() > 1) {
      MultiLine ll;
      for (auto& nf : getMainDirs()) {
        ll.push_back(nf.geom.getLine());
      }
      Polygon env = pbutil::geo::shrink(pbutil::geo::getOrientedEnvelopeAvg(ll), cd / 2).getPolygon();

      double incr = (bgeo::area(env) / bgeo::area(hull)) - 1;
      if (ll.size() < 5 || incr < 0.5) {
        hull = env;
      }
    }

    bgeo::buffer(hull, ret, distanceStrat, sideStrat, joinStrat, endStrat,
                 circleStrat);
  } else {
    // for two main dirs, take average
    std::vector<const PolyLine*> pols;

    PolyLine a = getMainDirs()[0].geom.getSegment(
        (cd / 2) / getMainDirs()[0].geom.getLength(),
        (getMainDirs()[0].geom.getLength() - cd / 2) /
            getMainDirs()[0].geom.getLength());

    PolyLine b = getMainDirs()[1].geom.getSegment(
        (cd / 2) / getMainDirs()[1].geom.getLength(),
        (getMainDirs()[1].geom.getLength() - cd / 2) /
            getMainDirs()[1].geom.getLength());

    assert(a.getLine().size() > 1);
    assert(b.getLine().size() > 1);

    if (dist(a.getLine()[0], b.getLine()[0]) >
        dist(a.getLine()[1], b.getLine()[0])) {
      a.reverse();
    }

    pols.push_back(&a);
    pols.push_back(&b);

    bgeo::buffer(PolyLine::average(pols).getLine(), ret, distanceStrat, sideStrat, joinStrat, endStrat,
                 circleStrat);
  }

  assert(ret.size() > 0);
  return ret[0];
}

// _____________________________________________________________________________
size_t Node::getNodeFrontPos(const NodeFront* a) const {
  for (size_t i = 0; i < _mainDirs.size(); ++i) {
    if (&_mainDirs[i] == a) return i;
  }

  return _mainDirs.size();
}

// _____________________________________________________________________________
void Node::addRouteConnException(const Route* r, const Edge* edgeA,
                                 const Edge* edgeB) {
  _routeConnExceptions[r][edgeA].insert(edgeB);
  // index the other direction also, will lead to faster lookups later on
  _routeConnExceptions[r][edgeB].insert(edgeA);
}

// _____________________________________________________________________________
bool Node::connOccurs(const Route* r, const Edge* edgeA,
                      const Edge* edgeB) const {
  const auto& i = _routeConnExceptions.find(r);
  if (i == _routeConnExceptions.end()) return true;
  const auto& ii = i->second.find(edgeA);
  if (ii == i->second.end()) return true;

  return ii->second.find(edgeB) == ii->second.end();
}

// _____________________________________________________________________________
Edge* Node::getEdge(const Node* other) const {
  for (auto e = _adjListOut.begin(); e != _adjListOut.end(); ++e) {
    Edge* eP = *e;

    if (eP->getTo() == other) {
      return eP;
    }
  }

  for (auto e = _adjListIn.begin(); e != _adjListIn.end(); ++e) {
    Edge* eP = *e;

    if (eP->getFrom() == other) {
      return eP;
    }
  }
}

// _____________________________________________________________________________
int Node::getCrossingPenalty(double inStatPen, double coef, bool adjpen) const {
  if (adjpen) coef *= _adjListIn.size() + _adjListOut.size();

  if (getStops().size() > 0) {
    return inStatPen * coef;
  }

  return coef;
}

// _____________________________________________________________________________
int Node::getSplittingPenalty(double inStatPen, double coef,
                              bool adjpen) const {
  if (adjpen) coef *= _adjListIn.size() + _adjListOut.size();

  if (getStops().size() > 0) {
    return inStatPen * coef;
  }

  return coef;
}
