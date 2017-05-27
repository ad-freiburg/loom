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
  RouteOccWithPos rop;

  rop = edge->getRouteOccWithPosUnder(r, c.find(edge)->second);

  if (rop.first) {
    return getTripPos(edge, rop.second, n == edge->getTo());
  }

  throw std::runtime_error("route does not occur in this edge");
}

// _____________________________________________________________________________
Point NodeFront::getTripPos(const Edge* e, size_t pos, bool inv) const {
  double p;
  if (!inv) {
    p = (e->getWidth() + e->getSpacing()) * pos + e->getWidth() / 2;
  } else {
    p = (e->getWidth() + e->getSpacing()) * (e->getCardinality() - 1 - pos) +
        e->getWidth() / 2;
  }

  // use interpolate here directly for speed
  return geom.interpolate(geom.getLine().front(), geom.getLine().back(), p);
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
double Node::getScore(const Configuration& c) const {
  std::vector<InnerGeometry> igs = getInnerGeometries(c, -1);
  size_t ret = 0;

  for (size_t i = 0; i < igs.size(); ++i) {
    for (size_t j = i + 1; j < igs.size(); ++j) {
      const InnerGeometry& iga = igs[i];
      const InnerGeometry& igb = igs[j];

      if (iga.from.front == igb.from.front && iga.slotFrom == igb.slotFrom) continue;
      if (iga.from.front == igb.to.front && iga.slotFrom == igb.slotTo) continue;
      if (iga.to.front == igb.to.front && iga.slotTo == igb.slotTo) continue;
      if (iga.to.front == igb.from.front && iga.slotTo == igb.slotFrom) continue;

      if (pbutil::geo::intersects(iga.geom.getLine().front(), iga.geom.getLine().back(), igb.geom.getLine().front(), igb.geom.getLine().back()) ||
          bgeo::distance(iga.geom.getLine(), igb.geom.getLine()) < 1) {
        ret++;
      }
    }
  }

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

        if (prec > 0) {
          ret.push_back(getInnerBezier(c, o, p, prec));
        } else {
          ret.push_back(getInnerStraightLine(c, o, p));
        }
      }

      processed[routeOcc.route].insert(&nf);
    }
  }

  return ret;
}

// _____________________________________________________________________________
InnerGeometry Node::getInnerStraightLine(
    const Configuration& c, const graph::Partner& partnerFrom,
    const graph::Partner& partnerTo) const {
  Point p = partnerFrom.front->getTripOccPos(partnerFrom.route, c);
  Point pp = partnerTo.front->getTripOccPos(partnerTo.route, c);

  size_t s = partnerFrom.edge->getRouteOccWithPosUnder(partnerFrom.route, c.find(partnerFrom.edge)->second).second;
  size_t ss = partnerTo.edge->getRouteOccWithPosUnder(partnerTo.route, c.find(partnerTo.edge)->second).second;

  return InnerGeometry(PolyLine(p, pp), partnerFrom, partnerTo, s, ss);
}

// _____________________________________________________________________________
InnerGeometry Node::getInnerBezier(const Configuration& cf,
                                   const graph::Partner& partnerFrom,
                                   const graph::Partner& partnerTo,
                                   double prec) const {
  InnerGeometry ret = getInnerStraightLine(cf, partnerFrom, partnerTo);
  Point p = ret.geom.getLine().front();
  Point pp = ret.geom.getLine().back();
  double d = pbutil::geo::dist(p, pp) / 2;

  Point b = p;
  Point c = pp;
  std::pair<double, double> slopeA, slopeB;

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

  b = Point(p.get<0>() + slopeA.first * d, p.get<1>() + slopeA.second * d);
  c = Point(pp.get<0>() + slopeB.first * d, pp.get<1>() + slopeB.second * d);

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
void Node::generateStationHull(double d) {
  _stationHull = getConvexFrontHull(d, true);
}

// _____________________________________________________________________________
Polygon Node::getStationHull() const {
  return _stationHull;
}

// _____________________________________________________________________________
Polygon Node::getConvexFrontHull(double d, bool rectangulize) const {
  MultiLine l;

  if (getMainDirs().size() != 2) {
    for (auto& nf : getMainDirs()) {
      // l.push_back(nf.geom.getLine());
      l.push_back(nf.geom.getSegment(
        (d / 2) /
            nf.geom.getLength(),
        (nf.geom.getLength() -
         d / 2) /
            nf.geom.getLength()).getLine());
    }
  } else {
    // for two main dirs, take average
    std::vector<const PolyLine*> pols;

    PolyLine a = getMainDirs()[0].geom.getSegment(
        (d / 2) /
            getMainDirs()[0].geom.getLength(),
        (getMainDirs()[0].geom.getLength() -
         d / 2) /
            getMainDirs()[0].geom.getLength());

    PolyLine b = getMainDirs()[1].geom.getSegment(
        (d / 2) /
            getMainDirs()[1].geom.getLength(),
        (getMainDirs()[1].geom.getLength() -
         d / 2) /
            getMainDirs()[1].geom.getLength());

    assert(a.getLine().size() > 1);
    assert(b.getLine().size() > 1);

    if (dist(a.getLine()[0], b.getLine()[0]) >
        dist(a.getLine()[1], b.getLine()[0])) {
      a.reverse();
    }

    pols.push_back(&a);
    pols.push_back(&b);
    l.push_back(PolyLine::average(pols).getLine());
  }

  MultiPolygon ret;
  double pointsPerCircle = 36;
  bgeo::strategy::buffer::distance_symmetric<double> distanceStrat(d);
  bgeo::strategy::buffer::join_round joinStrat(pointsPerCircle);
  bgeo::strategy::buffer::end_round endStrat(pointsPerCircle);
  bgeo::strategy::buffer::point_circle circleStrat(pointsPerCircle);
  bgeo::strategy::buffer::side_straight sideStrat;

  if (l.size() > 1) {
    Polygon hull;
    bgeo::convex_hull(l, hull);
    if (rectangulize && getMaxNodeFrontCardinality() > 1) {
      Polygon env = pbutil::geo::getOrientedEnvelopeAvg(l).getPolygon();
      double incr = (bgeo::area(env) / bgeo::area(hull)) - 1;
      if (bgeo::area(env) < d*d*36 && (l.size() < 5 || incr < 0.5)) {
        hull = env;
      }
    }

    bgeo::buffer(hull, ret, distanceStrat, sideStrat, joinStrat, endStrat,
                 circleStrat);
  } else {
    bgeo::buffer(l, ret, distanceStrat, sideStrat, joinStrat, endStrat,
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
