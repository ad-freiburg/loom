// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <set>
#include <string>
#include "json/json.hpp"
#include "shared/linegraph/Route.h"
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/TransitGraph.h"
#include "util/Misc.h"
#include "util/geo/BezierCurve.h"
#include "util/geo/Geo.h"

// we still need boost geometry for the polygon buffering, which is no
// yet implemented in util/Geo.h
#include <boost/geometry.hpp>
#include <boost/geometry/strategies/transform/matrix_transformers.hpp>

namespace bgeo = boost::geometry;

using util::geo::Point;
using util::geo::Box;
using util::geo::dist;
using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using shared::linegraph::Route;
using transitmapper::graph::OrderingConfig;
using transitmapper::graph::InnerGeometry;
using util::geo::BezierCurve;

// _____________________________________________________________________________
TransitGraph::TransitGraph() { _bbox = util::geo::minbox<double>(); }

// _____________________________________________________________________________
TransitGraph::~TransitGraph() {
  for (auto n : _nodes) {
    delete n;
  }
}

// _____________________________________________________________________________
const OrderingConfig& TransitGraph::getConfig() const { return _config; }

// _____________________________________________________________________________
void TransitGraph::setConfig(const OrderingConfig& c) { _config = c; }

// _____________________________________________________________________________
void TransitGraph::addNd(Node* n) {
  _nodes.insert(n);
  // expand the bounding box to hold this new node
  expandBBox(n->getPos());
}

// _____________________________________________________________________________
void TransitGraph::expandBBox(const DPoint& p) {
  _bbox = util::geo::extendBox(p, _bbox);
}

// _____________________________________________________________________________
Node* TransitGraph::getNodeById(const std::string& id) const {
  for (auto n : _nodes) {
    if (n->getId() == id) return n;
  }

  return 0;
}

// _____________________________________________________________________________
Edge* TransitGraph::addEdg(Node* from, Node* to, PolyLine<double> pl, double w,
                           double s) {
  if (from == to) return 0;
  Edge* e = getEdg(from, to);
  if (!e) {
    e = new Edge(from, to, pl, w, s);
    from->addEdg(e);
    to->addEdg(e);
    _bbox =
        util::geo::extendBox(util::geo::getBoundingBox(pl.getLine()), _bbox);
  }
  return e;
}

// _____________________________________________________________________________
void TransitGraph::delEdg(Node* from, Node* to) {
  Edge* toDel = getEdg(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdg(from, to));

  delete toDel;
}

// _____________________________________________________________________________
void TransitGraph::addRoute(const Route* r) {
  if (!getRoute(r->getId())) {
    _routes[r->getId()] = r;
  }
}

// _____________________________________________________________________________
const Route* TransitGraph::getRoute(const std::string& id) const {
  auto f = _routes.find(id);
  if (f == _routes.end()) return 0;

  return f->second;
}

// _____________________________________________________________________________
Edge* TransitGraph::getEdg(Node* from, Node* to) {
  for (auto e : from->getAdjListOut()) {
    if (e->getTo() == to) return e;
  }

  // also search in the opposite direction, we are handling an undirected
  // graph here
  for (auto e : from->getAdjListIn()) {
    if (e->getFrom() == to) return e;
  }

  return 0;
}

// _____________________________________________________________________________
const std::set<Node*>& TransitGraph::getNds() const { return _nodes; }

// _____________________________________________________________________________
std::set<Node*>* TransitGraph::getNds() { return &_nodes; }

// _____________________________________________________________________________
const DBox& TransitGraph::getBBox() const { return _bbox; }

// _____________________________________________________________________________
size_t TransitGraph::getNumNodes() const {
  return getNumNodes(true) + getNumNodes(false);
}

// _____________________________________________________________________________
size_t TransitGraph::getNumRoutes() const { return _routes.size(); }

// _____________________________________________________________________________
size_t TransitGraph::getMaxCardinality() const {
  size_t ret = 0;
  for (auto n : getNds()) {
    for (auto e : n->getAdjListOut()) {
      if (e->getCardinality() > ret) ret = e->getCardinality();
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::maxDeg() const {
  size_t ret = 0;
  for (auto n : getNds()) {
    if (n->getAdjListOut().size() + n->getAdjListIn().size() > ret) {
      ret = n->getAdjListOut().size() + n->getAdjListIn().size();
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumEdges() const {
  size_t ret = 0;

  for (auto n : getNds()) {
    ret += n->getAdjListOut().size();
  }

  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumNodes(bool topo) const {
  size_t ret = 0;
  for (auto n : _nodes) {
    if (n->getAdjListIn().size() + n->getAdjListOut().size() == 0) continue;
    if ((n->getStops().size() == 0) ^ !topo) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
bool TransitGraph::readFromJson(std::istream* s) {
  nlohmann::json j;
  (*s) >> j;

  if (j["type"] == "FeatureCollection") {
    // first pass, nodes
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = util::toString(props["id"]);

        // check if node already exists
        if (getNodeById(id)) continue;

        std::vector<double> coords = geom["coordinates"];

        Node* n = new Node(id, coords[0], coords[1]);

        shared::linegraph::Station i("", "", n->getPos());
        if (!props["station_id"].is_null() ||
            !props["station_label"].is_null()) {
          if (!props["station_id"].is_null())
            i.id = util::toString(props["station_id"]);
          if (!props["station_label"].is_null())
            i.name = util::toString(props["station_label"]);
          n->addStop(i);
        }

        addNd(n);
      }
    }

    // second pass, edges
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "LineString") {
        if (props["lines"].is_null() || props["lines"].size() == 0) continue;
        std::string from = util::toString(props["from"]);
        std::string to = util::toString(props["to"]);

        std::vector<std::vector<double>> coords = geom["coordinates"];

        PolyLine<double> pl;
        for (auto coord : coords) {
          double x = coord[0], y = coord[1];
          DPoint p(x, y);
          pl << p;
          expandBBox(p);
        }

        // TODO
        // pl.applyChaikinSmooth(_cfg->inputSmoothing);

        Node* fromN = getNodeById(from);
        Node* toN = getNodeById(to);

        if (!fromN) {
          continue;
        }

        if (!toN) {
          continue;
        }

        // TODO
        Edge* e = addEdg(fromN, toN, pl, 0, 0);

        assert(e);
        assert(getNodeById(from));
        assert(getNodeById(to));

        for (auto route : props["lines"]) {
          std::string id;
          if (!route["id"].is_null()) {
            id = util::toString(route["id"]);
          } else if (!route["label"].is_null()) {
            id = util::toString(route["label"]);
          } else if (!route["color"].is_null()) {
            id = route["color"];
          } else
            continue;

          const Route* r = getRoute(id);
          if (!r) {
            std::string label = route["label"].is_null() ? "" : route["label"];
            std::string color = route["color"];
            r = new Route(id, label, color);
            addRoute(r);
          }

          Node* dir = 0;

          if (!route["direction"].is_null()) {
            dir = getNodeById(util::toString(route["direction"]));
          }

          if (!route["style"].is_null()) {
            style::LineStyle ls;
            auto style = route["style"];
            std::string dashArray;
            if (!style["dash-array"].is_null()) {
              dashArray = style["dash-array"];
            }

            if (!style["css"].is_null()) {
              ls.setCss(style["css"]);
            }

            ls.setDashArray(dashArray);

            e->addRoute(r, dir, ls);
          } else {
            e->addRoute(r, dir);
          }
        }
      }
    }

    // third pass, exceptions (TODO: do this in the first part, store in some
    // data strcuture,
    //  add here!)
    for (auto feature : j["features"]) {
      auto props = feature["properties"];
      auto geom = feature["geometry"];
      if (geom["type"] == "Point") {
        std::string id = util::toString(props["id"]);

        Node* n = getNodeById(id);

        if (!n) continue;

        if (!props["excluded_line_conns"].is_null()) {
          for (auto excl : props["excluded_line_conns"]) {
            std::string rid = util::toString(excl["route"]);
            std::string nid1 = util::toString(excl["edge1_node"]);
            std::string nid2 = util::toString(excl["edge2_node"]);

            const Route* r = getRoute(rid);

            if (!r) {
              continue;
            }

            Node* n1 = getNodeById(nid1);
            Node* n2 = getNodeById(nid2);

            if (!n1) {
              continue;
            }

            if (!n2) {
              continue;
            }

            Edge* a = n->getEdg(n1);
            Edge* b = n->getEdg(n2);

            if (!a) {
              continue;
            }

            if (!b) {
              continue;
            }

            n->addRouteConnException(r, a, b);
          }
        }
      }
    }

  } else {
    return false;
  }

  return true;
}

// _____________________________________________________________________________
std::vector<InnerGeometry> TransitGraph::getInnerGeometries(
    const Node* n, const OrderingConfig& c, double prec) const {
  std::vector<InnerGeometry> ret;
  std::map<const Route*, std::set<const NodeFront*>> processed;

  for (size_t i = 0; i < n->getMainDirs().size(); ++i) {
    const NodeFront& nf = n->getMainDirs()[i];

    if (!c.count(nf.edge)) {
      std::cout << "No ordering for edge " << nf.edge << " found!" << std::endl;
      assert(false);
    }
    const std::vector<size_t>* ordering = &c.find(nf.edge)->second;

    for (size_t j : *ordering) {
      const RouteOccurance& routeOcc = (*nf.edge->getRoutes())[j];
      Partner o(&nf, nf.edge, routeOcc.route);

      std::vector<Partner> partners = n->getPartners(&nf, routeOcc);

      for (const Partner& p : partners) {
        if (processed[routeOcc.route].find(p.front) !=
            processed[routeOcc.route].end()) {
          continue;
        }

        auto is = getInnerStraightLine(n, c, o, p);
        if (is.geom.getLength() == 0) continue;

        if (prec > 0) {
          ret.push_back(getInnerBezier(n, c, o, p, prec));
        } else {
          ret.push_back(is);
        }
      }

      // handle lines where this node is the terminus
      if (partners.size() == 0) {
        auto is = getTerminusStraightLine(n, c, o);
        if (is.geom.getLength() > 0) {
          if (prec > 0) {
            ret.push_back(getTerminusBezier(n, c, o, prec));
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
InnerGeometry TransitGraph::getInnerBezier(const Node* n,
                                           const OrderingConfig& cf,
                                           const Partner& partnerFrom,
                                           const Partner& partnerTo,
                                           double prec) const {
  double EPSI = 0.001;
  InnerGeometry ret = getInnerStraightLine(n, cf, partnerFrom, partnerTo);
  DPoint p = ret.geom.getLine().front();
  DPoint pp = ret.geom.getLine().back();
  double d = util::geo::dist(p, pp);

  DPoint b = p;
  DPoint c = pp;
  std::pair<double, double> slopeA, slopeB;

  if (partnerFrom.front->edge->getGeom().getLength() <= 5) return ret;
  if (partnerTo.front->edge->getGeom().getLength() <= 5) return ret;

  if (d <= 5) return ret;

  if (partnerFrom.front->edge->getTo() == n) {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(
        partnerFrom.front->edge->getGeom().getLength() - 5,
        partnerFrom.front->edge->getGeom().getLength());
  } else {
    slopeA = partnerFrom.front->edge->getGeom().getSlopeBetweenDists(5, 0);
  }

  if (partnerTo.front->edge->getTo() == n) {
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

      da = (dar * 2 * (1 - ang) + avg * 2 * ang) / 2;
      db = (dbr * 2 * (1 - ang) + avg * 2 * ang) / 2;
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
InnerGeometry TransitGraph::getTerminusStraightLine(
    const Node* n, const OrderingConfig& c, const Partner& partnerFrom) const {
  DPoint p = partnerFrom.front->getTripOccPos(partnerFrom.route, c, false);
  DPoint pp = partnerFrom.front->getTripOccPos(partnerFrom.route, c, true);

  size_t s = partnerFrom.edge->getRoutePosUnder(
      partnerFrom.route, c.find(partnerFrom.edge)->second);
  size_t ss = partnerFrom.edge->getRoutePosUnder(
      partnerFrom.route, c.find(partnerFrom.edge)->second);

  return InnerGeometry(PolyLine<double>(p, pp), partnerFrom, Partner(), s, ss);
}

// _____________________________________________________________________________
InnerGeometry TransitGraph::getInnerStraightLine(
    const Node* n, const OrderingConfig& c, const Partner& partnerFrom,
    const Partner& partnerTo) const {
  DPoint p = partnerFrom.front->getTripOccPos(partnerFrom.route, c, false);
  DPoint pp = partnerTo.front->getTripOccPos(partnerTo.route, c, false);

  size_t s = partnerFrom.edge->getRoutePosUnder(
      partnerFrom.route, c.find(partnerFrom.edge)->second);
  size_t ss = partnerTo.edge->getRoutePosUnder(partnerTo.route,
                                               c.find(partnerTo.edge)->second);

  return InnerGeometry(PolyLine<double>(p, pp), partnerFrom, partnerTo, s, ss);
}

// _____________________________________________________________________________
InnerGeometry TransitGraph::getTerminusBezier(const Node* n,
                                              const OrderingConfig& cf,
                                              const Partner& partnerFrom,
                                              double prec) const {
  InnerGeometry ret = getTerminusStraightLine(n, cf, partnerFrom);
  DPoint p = ret.geom.getLine().front();
  DPoint pp = ret.geom.getLine().back();
  double d = util::geo::dist(p, pp) / 2;

  DPoint b = p;
  DPoint c = pp;
  std::pair<double, double> slopeA;

  assert(partnerFrom.front->edge->getGeom().getLength() > 5);

  if (partnerFrom.front->edge->getTo() == n) {
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
Polygon<double> TransitGraph::getConvexFrontHull(
    const Node* n, double d, bool rectangulize,
    bool simpleRenderForTwoEdgeNodes) const {
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

  if (!simpleRenderForTwoEdgeNodes || n->getMainDirs().size() != 2) {
    MultiLine<double> l;
    for (auto& nf : n->getMainDirs()) {
      l.push_back(
          nf.geom
              .getSegment((cd / 2) / nf.geom.getLength(),
                          (nf.geom.getLength() - cd / 2) / nf.geom.getLength())
              .getLine());
    }

    Polygon<double> hull = util::geo::convexHull(l);

    if (rectangulize && n->getMaxNodeFrontCardinality() > 1) {
      MultiLine<double> ll;
      for (auto& nf : n->getMainDirs()) {
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

    PolyLine<double> a = n->getMainDirs()[0].geom.getSegment(
        (cd / 2) / n->getMainDirs()[0].geom.getLength(),
        (n->getMainDirs()[0].geom.getLength() - cd / 2) /
            n->getMainDirs()[0].geom.getLength());

    PolyLine<double> b = n->getMainDirs()[1].geom.getSegment(
        (cd / 2) / n->getMainDirs()[1].geom.getLength(),
        (n->getMainDirs()[1].geom.getLength() - cd / 2) /
            n->getMainDirs()[1].geom.getLength());

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
Polygon<double> TransitGraph::getStationHull(const Node* n, double d,
                                             bool simple) const {
  if (n->getMainDirs().size() == 0) return Polygon<double>();
  return getConvexFrontHull(n, d, true, simple);
}
