// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <algorithm>
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/Node.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/Geo.h"

// we still need boost geometry for the polygon buffering, which is no
// yet implemented in util/Geo.h
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>

using namespace transitmapper;
using namespace graph;

namespace bgeo = boost::geometry;

using util::geo::Point;
using util::geo::Polygon;
using util::geo::BezierCurve;
using util::geo::PolyLine;
using util::geo::Line;

// _____________________________________________________________________________
DPoint NodeFront::getTripOccPos(const Route* r, const OrderingConfig& c) const {
  return getTripOccPos(r, c, false);
}

// _____________________________________________________________________________
DPoint NodeFront::getTripOccPos(const Route* r, const OrderingConfig& c,
                                bool origGeom) const {
  RouteOccWithPos rop;

  assert(c.find(edge) != c.end());

  rop = edge->getRouteWithPosUnder(r, c.find(edge)->second);

  if (rop.first) {
    return getTripPos(edge, rop.second, n == edge->getTo(), origGeom);
  }

  throw std::runtime_error("route does not occur in this edge");
}

// _____________________________________________________________________________
DPoint NodeFront::getTripPos(const Edge* e, size_t pos, bool inv) const {
  return getTripPos(e, pos, inv, false);
}

// _____________________________________________________________________________
double NodeFront::getOutAngle() const {
  double checkDist = 10;
  if (edge->getFrom() == n) {
    return angBetween(n->getPos(), edge->getGeom().getPointAtDist(checkDist).p);
  } else {
    return angBetween(
        n->getPos(),
        edge->getGeom()
            .getPointAtDist(edge->getGeom().getLength() - checkDist)
            .p);
  }
}

// _____________________________________________________________________________
DPoint NodeFront::getTripPos(const Edge* e, size_t pos, bool inv,
                             bool origG) const {
  double p;
  if (!inv) {
    p = (e->getWidth() + e->getSpacing()) * pos + e->getWidth() / 2;
  } else {
    p = (e->getWidth() + e->getSpacing()) * (e->getCardinality() - 1 - pos) +
        e->getWidth() / 2;
  }
  // use interpolate here directly for speed
  if (origG) {
    return origGeom.interpolate(origGeom.getLine().front(),
                                origGeom.getLine().back(), p);
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
Node::Node(const std::string& id, DPoint pos) : _id(id), _pos(pos) {}

// _____________________________________________________________________________
Node::Node(const std::string& id, double x, double y) : _id(id), _pos(x, y) {}

// _____________________________________________________________________________
Node::Node(const std::string& id, DPoint pos, StationInfo s)
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
const DPoint& Node::getPos() const { return _pos; }

// _____________________________________________________________________________
void Node::setPos(const DPoint& p) { _pos = p; }

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
         nf.edge->getCtdRoutesIn(this, ro.route, ro.direction, f->edge)) {
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
std::vector<InnerGeometry> Node::getInnerGeometries(const OrderingConfig& c,
                                                    double prec) const {
  std::vector<InnerGeometry> ret;
  std::map<const Route*, std::set<const NodeFront*>> processed;

  for (size_t i = 0; i < getMainDirs().size(); ++i) {
    const NodeFront& nf = getMainDirs()[i];

    if (!c.count(nf.edge)) {
      std::cout << "No ordering for edge " << nf.edge << " found!" << std::endl;
      assert(false);
    }
    const std::vector<size_t>* ordering = &c.find(nf.edge)->second;

    for (size_t j : *ordering) {
      const RouteOccurance& routeOcc = (*nf.edge->getRoutes())[j];
      Partner o(&nf, nf.edge, routeOcc.route);

      std::vector<Partner> partners = getPartners(&nf, routeOcc);

      for (const Partner& p : partners) {
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
InnerGeometry Node::getTerminusStraightLine(const OrderingConfig& c,
                                            const Partner& partnerFrom) const {
  DPoint p = partnerFrom.front->getTripOccPos(partnerFrom.route, c, false);
  DPoint pp = partnerFrom.front->getTripOccPos(partnerFrom.route, c, true);

  size_t s = partnerFrom.edge
                 ->getRouteWithPosUnder(partnerFrom.route,
                                        c.find(partnerFrom.edge)->second)
                 .second;
  size_t ss = partnerFrom.edge
                  ->getRouteWithPosUnder(partnerFrom.route,
                                         c.find(partnerFrom.edge)->second)
                  .second;

  return InnerGeometry(PolyLine<double>(p, pp), partnerFrom, Partner(), s, ss);
}

// _____________________________________________________________________________
InnerGeometry Node::getInnerStraightLine(const OrderingConfig& c,
                                         const Partner& partnerFrom,
                                         const Partner& partnerTo) const {
  DPoint p = partnerFrom.front->getTripOccPos(partnerFrom.route, c);
  DPoint pp = partnerTo.front->getTripOccPos(partnerTo.route, c);

  size_t s = partnerFrom.edge
                 ->getRouteWithPosUnder(partnerFrom.route,
                                        c.find(partnerFrom.edge)->second)
                 .second;
  size_t ss = partnerTo.edge
                  ->getRouteWithPosUnder(partnerTo.route,
                                         c.find(partnerTo.edge)->second)
                  .second;

  return InnerGeometry(PolyLine<double>(p, pp), partnerFrom, partnerTo, s, ss);
}

// _____________________________________________________________________________
InnerGeometry Node::getTerminusBezier(const OrderingConfig& cf,
                                      const Partner& partnerFrom,
                                      double prec) const {
  InnerGeometry ret = getTerminusStraightLine(cf, partnerFrom);
  DPoint p = ret.geom.getLine().front();
  DPoint pp = ret.geom.getLine().back();
  double d = util::geo::dist(p, pp) / 2;

  DPoint b = p;
  DPoint c = pp;
  std::pair<double, double> slopeA;

  assert(partnerFrom.front->edge->getGeom().getLength() > 5);

  if (partnerFrom.front->edge->getTo() == this) {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(
        partnerFrom.front->edge->getGeom().getLength() - 5,
        partnerFrom.front->edge->getGeom().getLength());
  } else {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(5, 0);
  }

  b = DPoint(p.getX() + slopeA.first * d, p.getY() + slopeA.second * d);
  c = DPoint(pp.getX(), pp.getY());

  BezierCurve<double> bc(p, b, c, pp);
  ret.geom = bc.render(prec);

  return ret;
}

// _____________________________________________________________________________
InnerGeometry Node::getInnerBezier(const OrderingConfig& cf,
                                   const Partner& partnerFrom,
                                   const Partner& partnerTo,
                                   double prec) const {
  double EPSI = 0.001;
  InnerGeometry ret = getInnerStraightLine(cf, partnerFrom, partnerTo);
  DPoint p = ret.geom.getLine().front();
  DPoint pp = ret.geom.getLine().back();
  double d = util::geo::dist(p, pp);

  DPoint b = p;
  DPoint c = pp;
  std::pair<double, double> slopeA, slopeB;

  if (partnerFrom.front->edge->getGeom().getLength() <= 5) return ret;
  if (partnerTo.front->edge->getGeom().getLength() <= 5) return ret;

  if (d <= 5) return ret;

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

  DPoint pa =
      DPoint(p.getX() - slopeA.first * d * 2, p.getY() - slopeA.second * d * 2);
  DPoint pb = DPoint(pp.getX() - slopeB.first * d * 2,
                     pp.getY() - slopeB.second * d * 2);
  DPoint ppa =
      DPoint(p.getX() + slopeA.first * d * 2, p.getY() + slopeA.second * d * 2);
  DPoint ppb = DPoint(pp.getX() + slopeB.first * d * 2,
                      pp.getY() + slopeB.second * d * 2);

  if (d > EPSI && util::geo::intersects(pa, ppa, pb, ppb)) {
    DPoint isect = util::geo::intersection(pa, ppa, pb, ppb);

    if (!std::isnan(isect.getX()) && !std::isnan(isect.getY())) {
      double degAng = util::geo::innerProd(isect, pa, ppb);
      double ang = cos(degAng / (180 / M_PI));

      ang = ang * ang;

      if (std::isnan(ang)) ang = 1;

      if (std::max(partnerFrom.edge->getCardinality(),
                   partnerTo.edge->getCardinality()) > 1) {
        double fac = fabs((double)((int)partnerFrom.edge->getCardinality() -
                                   (int)partnerTo.edge->getCardinality())) /
                     (double)(std::max(partnerFrom.edge->getCardinality(),
                                       partnerTo.edge->getCardinality()) -
                              1);

        fac = pow(fac, 2);
        ang = pow(ang, fac);
      }

      double dar = util::geo::dist(isect, p);
      double dbr = util::geo::dist(isect, pp);

      double avg = (dar + dbr) / 2;

      da = (dar * 2 * (1 - ang) + avg * 2 * (ang)) / 2;
      db = (dbr * 2 * (1 - ang) + avg * 2 * (ang)) / 2;
    }
  }

  if ((da + db) < EPSI || (da + db) > d * 2) {
    da = 1;
    db = 1;
  }

  b = DPoint(p.getX() + slopeA.first * d * (da / (da + db)),
             p.getY() + slopeA.second * d * (da / (da + db)));
  c = DPoint(pp.getX() + slopeB.first * d * (db / (da + db)),
             pp.getY() + slopeB.second * d * (db / (da + db)));

  BezierCurve<double> bc(p, b, c, pp);
  ret.geom = bc.render(prec);

  return ret;
}

// _____________________________________________________________________________
size_t Node::getConnCardinality() const {
  size_t ret = 0;
  std::map<const Route*, std::set<const NodeFront*>> processed;

  for (size_t i = 0; i < getMainDirs().size(); ++i) {
    const NodeFront& nf = getMainDirs()[i];

    for (size_t j = 0; j < nf.edge->getCardinality(); j++) {
      const RouteOccurance& routeOcc = (*nf.edge->getRoutes())[j];

      std::vector<Partner> partners = getPartners(&nf, routeOcc);

      for (const Partner& p : partners) {
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
Polygon<double> Node::getStationHull() const { return _stationHull; }

// _____________________________________________________________________________
Polygon<double> Node::getConvexFrontHull(
    double d, bool rectangulize, bool simpleRenderForTwoEdgeNodes) const {
  double cd = d;

  typedef bgeo::model::point<double, 2, bgeo::cs::cartesian> BoostPoint;
  typedef bgeo::model::linestring<BoostPoint> BoostLine;
  typedef bgeo::model::polygon<BoostPoint> BoostPoly;
  typedef bgeo::model::multi_polygon<BoostPoly> BoostMultiPoly;

  BoostMultiPoly ret;
  double pointsPerCircle = 36;
  bgeo::strategy::buffer::distance_symmetric<double> distanceStrat(d);
  bgeo::strategy::buffer::join_round joinStrat(pointsPerCircle);
  bgeo::strategy::buffer::end_round endStrat(pointsPerCircle);
  bgeo::strategy::buffer::point_circle circleStrat(pointsPerCircle);
  bgeo::strategy::buffer::side_straight sideStrat;

  if (!simpleRenderForTwoEdgeNodes || getMainDirs().size() != 2) {
    MultiLine<double> l;
    for (auto& nf : getMainDirs()) {
      l.push_back(
          nf.geom
              .getSegment((cd / 2) / nf.geom.getLength(),
                          (nf.geom.getLength() - cd / 2) / nf.geom.getLength())
              .getLine());
    }

    Polygon<double> hull = util::geo::convexHull(l);

    if (rectangulize && getMaxNodeFrontCardinality() > 1) {
      MultiLine<double> ll;
      for (auto& nf : getMainDirs()) {
        ll.push_back(nf.geom.getLine());
      }
      Polygon<double> env = util::geo::convexHull(
          util::geo::shrink(util::geo::getOrientedEnvelopeAvg(ll), cd / 2));

      double incr = (util::geo::area(env) / util::geo::area(hull)) - 1;
      if (ll.size() < 5 || incr < 0.5) {
        hull = env;
      }
    }

    BoostPoly hullBgeo;
    for (const auto& p : hull.getOuter())
      hullBgeo.outer().push_back({p.getX(), p.getY()});
    hullBgeo.outer().push_back(
        {hull.getOuter().front().getX(), hull.getOuter().front().getY()});

    // boost geometry expects polygons in clockwise fashion
    bgeo::correct(hullBgeo);

    bgeo::buffer(hullBgeo, ret, distanceStrat, sideStrat, joinStrat, endStrat,
                 circleStrat);

  } else {
    // for two main dirs, take average
    std::vector<const PolyLine<double>*> pols;

    PolyLine<double> a = getMainDirs()[0].geom.getSegment(
        (cd / 2) / getMainDirs()[0].geom.getLength(),
        (getMainDirs()[0].geom.getLength() - cd / 2) /
            getMainDirs()[0].geom.getLength());

    PolyLine<double> b = getMainDirs()[1].geom.getSegment(
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

    auto avg = PolyLine<double>::average(pols).getLine();

    BoostLine lineBgeo;
    for (const auto& p : avg) {
      lineBgeo.push_back({p.getX(), p.getY()});
    }

    bgeo::buffer(lineBgeo, ret, distanceStrat, sideStrat, joinStrat, endStrat,
                 circleStrat);
  }

  assert(ret.size() > 0);

  Polygon<double> retPoly;
  for (const auto& p : ret[0].outer()) {
    retPoly.getOuter().push_back({p.get<0>(), p.get<1>()});
  }

  return retPoly;
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

  return 0;
}
