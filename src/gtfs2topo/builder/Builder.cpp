// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include "./Builder.h"
#include "pbutil/geo/Geo.h"
#include "pbutil/log/Log.h"

using namespace gtfs2topo;
using namespace graph;
using namespace gtfsparser;
using namespace gtfs;

using pbutil::geo::Point;
using pbutil::geo::PolyLine;
using pbutil::geo::SharedSegments;

// _____________________________________________________________________________
Builder::Builder(const config::Config* cfg) : _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
}

// _____________________________________________________________________________
void Builder::consume(const Feed& f, Graph* g) {
  // TODO: make this stuff configurable

  uint8_t AGGREGATE_STOPS = _cfg->stationAggrLevel;

  size_t cur = 0;
  for (auto t = f.tripsBegin(); t != f.tripsEnd(); ++t) {
    cur++;
    if (t->second->getStopTimes().size() < 2) continue;

    if (t->second->getRoute()->getType() == gtfs::Route::TYPE::BUS) continue;
    // if (t->second->getRoute()->getType() != gtfs::Route::TYPE::RAIL) continue;
    // if (t->second->getRoute()->getType() != gtfs::Route::TYPE::TRAM) continue;
    if (!checkTripSanity(t->second)) continue;

    auto st = t->second->getStopTimes().begin();

    StopTime prev = *st;
    const Edge* prevEdge = 0;
    addStop(prev.getStop(), AGGREGATE_STOPS, g);
    ++st;

    for (; st != t->second->getStopTimes().end(); ++st) {
      const StopTime& cur = *st;

      Node* fromNode = g->getNodeByStop(prev.getStop(), AGGREGATE_STOPS);

      Node* toNode = addStop(cur.getStop(), AGGREGATE_STOPS, g);

      if (fromNode == toNode) continue;

      Edge* exE = g->getEdge(fromNode, toNode);

      if (!exE) {
        exE = g->addEdge(fromNode, toNode);
      }

      Node* directionNode = toNode;

      if (!exE->addTrip(t->second, directionNode)) {
        std::pair<bool, PolyLine> edgeGeom;
        if (AGGREGATE_STOPS) {
          Stop* frs = prev.getStop()->getParentStation()
                          ? prev.getStop()->getParentStation()
                          : prev.getStop();
          Stop* tos = cur.getStop()->getParentStation()
                          ? cur.getStop()->getParentStation()
                          : cur.getStop();
          edgeGeom = getSubPolyLine(
              frs, tos, t->second, prev.getShapeDistanceTravelled(),
              cur.getShapeDistanceTravelled(), g->getProjection());
        } else {
          edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second,
                                    prev.getShapeDistanceTravelled(),
                                    cur.getShapeDistanceTravelled(),
                                    g->getProjection());
        }

        // only take geometries that could be found using the
        // shape geometry, not some fallback
        if (edgeGeom.first) {
          if (prevEdge) {
            fromNode->connOccurs(t->second->getRoute(), prevEdge, exE);
          }

          exE->addTrip(t->second, edgeGeom.second, directionNode);
        }
      }
      prev = cur;
      prevEdge = exE;
    }
  }
}

// _____________________________________________________________________________
Point Builder::getProjectedPoint(double lat, double lng, projPJ p) const {
  double x = lng;
  double y = lat;
  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;

  pj_transform(_mercProj, p, 1, 1, &x, &y, 0);

  return Point(x, y);
}

// _____________________________________________________________________________
void Builder::simplify(Graph* g) {
  // try to merge both-direction edges into a single one

  for (auto n : *g->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      e->simplify();
    }
  }
}

// _____________________________________________________________________________
ShrdSegWrap Builder::getNextSharedSegment(Graph* g, bool final) const {
  int i = 0;
  double dmin = _cfg->maxAggrDistance;
  for (auto n : *g->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      i++;
      if (_indEdges.find(e) != _indEdges.end() ||
          e->getEdgeTripGeoms()->size() != 1) {
        continue;
      }
      // TODO: outfactor this _______
      for (auto nt : *g->getNodes()) {
        for (auto toTest : nt->getAdjListOut()) {
          // TODO: only check edges with a SINGLE geometry atm, see also check
          // above
          if (_indEdges.find(toTest) != _indEdges.end() ||
              toTest->getEdgeTripGeoms()->size() != 1) {
            continue;
          }

          if (e != toTest) {
            double dmax = 5;

            if (final) {
              // set dmax according to the width of the etg
              dmax = fmax(
                  dmin, e->getEdgeTripGeoms()->front().getCardinality() *
                              (dmin / 2) +
                          toTest->getEdgeTripGeoms()->front().getCardinality() *
                              (dmin / 2));
            }

            SharedSegments s =
                e->getEdgeTripGeoms()->front().getGeom().getSharedSegments(
                    toTest->getEdgeTripGeoms()->front().getGeom(), dmax);

            if (s.segments.size() > 0) {
              _pEdges[e]++;
              if (_pEdges[e] > 20) {
                _indEdges.insert(e);
              }
              // LOG(DEBUG) << _indEdges.size() << " / " << i << std::endl;
              return ShrdSegWrap(e, toTest, s.segments.front());
            }
          }
        }
      }

      // nothing overlapping found for this edge anymore, mark as processed
      _indEdges.insert(e);
      // _________
    }
  }

  return ShrdSegWrap();
}

// _____________________________________________________________________________
bool Builder::lineCrossesAtNode(const Node* a, const Edge* ea,
                                const Edge* eb) const {
  if (!(ea->getTo() == a || ea->getFrom() == a) ||
      !(eb->getTo() == a || eb->getFrom() == a))
    return false;

  double lookDist = 100;
  double minInAngle = 0.2;

  const EdgeTripGeom& e = ea->getEdgeTripGeoms().front();
  const EdgeTripGeom& f = eb->getEdgeTripGeoms().front();

  double angleE = 0;
  double angleF = 0;
  double angleEPartner = 0;
  double angleFPartner = 0;

  if (e.getGeomDir() == a) {
    // look at end of geom
    Point a = e.getGeom().getLine().back();
    Point b = e.getGeom().getPointAtDist(e.getGeom().getLength() - lookDist).p;
    angleE = pbutil::geo::angBetween(a, b);
  } else {
    // look at front of geom
    Point a = e.getGeom().getLine().front();
    Point b = e.getGeom().getPointAtDist(lookDist).p;
    angleE = pbutil::geo::angBetween(a, b);
  }

  if (f.getGeomDir() == a) {
    // look at end of geom
    Point a = f.getGeom().getLine().back();
    Point b = f.getGeom().getPointAtDist(f.getGeom().getLength() - lookDist).p;
    angleF = pbutil::geo::angBetween(a, b);
  } else {
    // look at front of geom
    Point a = f.getGeom().getLine().front();
    Point b = f.getGeom().getPointAtDist(lookDist).p;
    angleF = pbutil::geo::angBetween(a, b);
  }

  double inDiff = fabs(angleE - angleF);

  if (inDiff < minInAngle) return false;

  bool eFound = false;
  bool fFound = false;

  // now lets get some edge that is not e or f and has all of e's routes

  for (auto edge : a->getAdjList()) {
    if (edge == ea || edge == eb) continue;
    if (edge->getEdgeTripGeoms()->size() == 0) continue;
    const EdgeTripGeom& etg = edge->getEdgeTripGeoms()->front();

    if (etg.routeEquivalent(e)) {
      if (etg.getGeomDir() != a) {
        // look at end of geom
        Point a = etg.getGeom().getLine().back();
        Point b = etg.getGeom()
                      .getPointAtDist(etg.getGeom().getLength() - lookDist)
                      .p;
        angleEPartner = pbutil::geo::angBetween(a, b);
      } else {
        // look at front of geom
        Point a = etg.getGeom().getLine().front();
        Point b = etg.getGeom().getPointAtDist(lookDist).p;
        angleEPartner = pbutil::geo::angBetween(a, b);
      }
      eFound = true;
      break;
    }
  }

  for (auto edge : a->getAdjList()) {
    if (edge == ea || edge == eb) continue;
    if (edge->getEdgeTripGeoms()->size() == 0) continue;
    const EdgeTripGeom& etg = edge->getEdgeTripGeoms()->front();

    if (etg.routeEquivalent(f)) {
      if (etg.getGeomDir() != a) {
        // look at end of geom
        Point a = etg.getGeom().getLine().back();
        Point b = etg.getGeom()
                      .getPointAtDist(etg.getGeom().getLength() - lookDist)
                      .p;
        angleFPartner = pbutil::geo::angBetween(a, b);
      } else {
        // look at front of geom
        Point a = etg.getGeom().getLine().front();
        Point b = etg.getGeom().getPointAtDist(lookDist).p;
        angleFPartner = pbutil::geo::angBetween(a, b);
      }
      fFound = true;
      break;
    }
  }

  if (!eFound || !fFound) return false;

  double outDiff = fabs(angleEPartner - angleFPartner);

  if (outDiff < minInAngle) return false;

  double aDiff = fabs(angleE - angleEPartner);
  double bDiff = fabs(angleF - angleFPartner);

  if (aDiff < 0.2 && bDiff < 0.2) {
    return true;
  }

  return false;
}

// _____________________________________________________________________________
bool Builder::createTopologicalNodes(Graph* g, bool final) {
  ShrdSegWrap w;
  _indEdges.clear();
  _pEdges.clear();
  bool found = false;

  while ((w = getNextSharedSegment(g, final)).e) {
    const EdgeTripGeom& curEdgeGeom = w.e->getEdgeTripGeoms()->front();
    const EdgeTripGeom& compEdgeGeom = w.f->getEdgeTripGeoms()->front();

    const LinePoint& eap = w.s.first.first;
    const LinePoint& ebp = w.s.second.first;
    const LinePoint& fap = w.s.first.second;
    const LinePoint& fbp = w.s.second.second;

    PolyLine ea = curEdgeGeom.getGeom().getSegment(0, eap.totalPos);
    PolyLine ec = curEdgeGeom.getGeom().getSegment(ebp.totalPos, 1);

    PolyLine fa;
    PolyLine fc;

    assert(w.f->getTo() == compEdgeGeom.getGeomDir());

    if (fap.totalPos > fbp.totalPos) {
      fa = compEdgeGeom.getGeom().getSegment(fap.totalPos, 1);
      fc = compEdgeGeom.getGeom().getSegment(0, fbp.totalPos);
    } else {
      fa = compEdgeGeom.getGeom().getSegment(0, fap.totalPos);
      fc = compEdgeGeom.getGeom().getSegment(fbp.totalPos, 1);
    }

    PolyLine ab = getAveragedFromSharedSeg(w);

    // new nodes at the start and end of the shared segment
    Node* a = 0;
    Node* b = 0;

    double maxSnapDist = 20;

    if (ea.getLength() < maxSnapDist) {
      a = w.e->getFrom();
    }

    if (ec.getLength() < maxSnapDist) {
      b = w.e->getTo();
    }

    if (fa.getLength() < maxSnapDist) {
      a = !(fap.totalPos > fbp.totalPos) ? w.f->getFrom() : w.f->getTo();
    }

    if (fc.getLength() < maxSnapDist) {
      b = !(fap.totalPos > fbp.totalPos) ? w.f->getTo() : w.f->getFrom();
    }

    if (!a) a = new Node(w.s.first.first.p);
    if (!b) b = new Node(w.s.second.first.p);

    if (a == b) {
      continue;
    }

    if (lineCrossesAtNode(a, w.e, w.f)) continue;
    if (lineCrossesAtNode(b, w.e, w.f)) continue;

    EdgeTripGeom eaEdgeGeom(ea, a);
    EdgeTripGeom abEdgeGeom(ab, b);
    EdgeTripGeom ecEdgeGeom(ec, w.e->getTo());

    const Node* faDir = 0;
    const Node* fcDir = 0;

    if (fap.totalPos > fbp.totalPos) {
      faDir = w.f->getTo();
      fcDir = b;
    } else {
      faDir = a;
      fcDir = w.f->getTo();
    }

    EdgeTripGeom faEdgeGeom(fa, faDir);
    EdgeTripGeom fcEdgeGeom(fc, fcDir);

    for (const TripOccurance& r : curEdgeGeom.getTripsUnordered()) {
      for (auto& t : r.trips) {
        if (!r.direction) {
          eaEdgeGeom.addTrip(t, 0);
          abEdgeGeom.addTrip(t, 0);
          ecEdgeGeom.addTrip(t, 0);
        } else if (r.direction == w.e->getTo()) {
          eaEdgeGeom.addTrip(t, a);
          abEdgeGeom.addTrip(t, b);
          ecEdgeGeom.addTrip(t, w.e->getTo());
        } else {
          eaEdgeGeom.addTrip(t, w.e->getFrom());
          abEdgeGeom.addTrip(t, a);
          ecEdgeGeom.addTrip(t, b);
        }
      }
    }

    for (const TripOccurance r : compEdgeGeom.getTripsUnordered()) {
      for (auto& t : r.trips) {
        if (!r.direction) {
          faEdgeGeom.addTrip(t, 0);
          abEdgeGeom.addTrip(t, 0);
          fcEdgeGeom.addTrip(t, 0);
        } else if ((r.direction == w.f->getTo() &&
                    !(fap.totalPos > fbp.totalPos)) ||
                   (r.direction != w.f->getTo() &&
                    (fap.totalPos > fbp.totalPos))) {
          faEdgeGeom.addTrip(t, a);
          abEdgeGeom.addTrip(t, b);
          fcEdgeGeom.addTrip(t, w.f->getTo());
        } else {
          faEdgeGeom.addTrip(t, w.f->getFrom());
          abEdgeGeom.addTrip(t, a);
          fcEdgeGeom.addTrip(t, b);
        }
      }
    }

    for (auto to : curEdgeGeom.getTripsUnordered()) {
      assert(abEdgeGeom.containsRoute(to.route));
    }

    for (auto to : compEdgeGeom.getTripsUnordered()) {
      assert(abEdgeGeom.containsRoute(to.route));
    }

    Node* wefrom = w.e->getFrom();
    Node* weto = w.e->getTo();
    Node* wffrom = w.f->getFrom();
    Node* wfto = w.f->getTo();

    // delete old edges
    g->deleteEdge(w.e->getFrom(), w.e->getTo());
    g->deleteEdge(w.f->getFrom(), w.f->getTo());

    // add new edges
    g->addNode(a);
    g->addNode(b);
    Edge* eaE = g->addEdge(wefrom, a);
    Edge* abE = g->addEdge(a, b);
    Edge* ebE = g->addEdge(b, weto);

    if (eaE) {
      wefrom->replaceEdgeInConnections(w.e, eaE);
    } else {
      assert(a == wefrom);
      a->replaceEdgeInConnections(w.e, abE);
    }

    if (ebE) {
      weto->replaceEdgeInConnections(w.e, ebE);
    } else {
      assert(b == weto);
      b->replaceEdgeInConnections(w.e, abE);
    }

    Edge* faE = 0;
    Edge* fbE = 0;

    if (fap.totalPos > fbp.totalPos) {
      faE = g->addEdge(a, wfto);
      fbE = g->addEdge(wffrom, b);

      if (faE) {
        wfto->replaceEdgeInConnections(w.f, faE);
      } else {
        assert(a == wfto);
        a->replaceEdgeInConnections(w.f, abE);
      }

      if (fbE) {
        wffrom->replaceEdgeInConnections(w.f, fbE);
      } else {
        assert(b == wffrom);
        b->replaceEdgeInConnections(w.f, abE);
      }
    } else {
      faE = g->addEdge(wffrom, a);
      fbE = g->addEdge(b, wfto);

      if (faE) {
        wffrom->replaceEdgeInConnections(w.f, faE);
      } else {
        assert(a == wffrom);
        a->replaceEdgeInConnections(w.f, abE);
      }

      if (fbE) {
        wfto->replaceEdgeInConnections(w.f, fbE);
      } else {
        assert(b == wfto);
        b->replaceEdgeInConnections(w.f, abE);
      }
    }

    if (abE) {
      abE->addEdgeTripGeom(abEdgeGeom);
      abE->simplify();
    } else {
      // we use abE below without checking!
      assert(false);
    }

    if (eaE) {
      eaE->addEdgeTripGeom(eaEdgeGeom);
      eaE->simplify();
      a->sewConnectionsTogether(eaE, abE);
    }
    if (ebE) {
      ebE->addEdgeTripGeom(ecEdgeGeom);
      ebE->simplify();
      b->sewConnectionsTogether(abE, ebE);
    }

    if (faE) {
      faE->addEdgeTripGeom(faEdgeGeom);
      faE->simplify();
      a->sewConnectionsTogether(faE, abE);
    }
    if (fbE) {
      fbE->addEdgeTripGeom(fcEdgeGeom);
      fbE->simplify();
      b->sewConnectionsTogether(abE, fbE);
    }

    found = true;
  }

  return found;
}

// _____________________________________________________________________________
std::pair<bool, PolyLine> Builder::getSubPolyLine(Stop* a, Stop* b, Trip* t,
                                                  double distA, double distB,
                                                  projPJ proj) {
  Point ap = getProjectedPoint(a->getLat(), a->getLng(), proj);
  Point bp = getProjectedPoint(b->getLat(), b->getLng(), proj);

  if (!t->getShape()) {
    return std::pair<bool, PolyLine>(false, PolyLine(ap, bp));
  }

  double totalTripDist = t->getShape()->getPoints().rbegin()->travelDist -
                         t->getShape()->getPoints().begin()->travelDist;

  auto pl = _polyLines.find(t->getShape());
  if (pl == _polyLines.end()) {
    // generate polyline for this shape
    pl = _polyLines
             .insert(
                 std::pair<gtfs::Shape*, PolyLine>(t->getShape(), PolyLine()))
             .first;

    for (const auto& sp : t->getShape()->getPoints()) {
      pl->second << getProjectedPoint(sp.lat, sp.lng, proj);
    }

    pl->second.simplify(20);
    pl->second.smoothenOutliers(50);
    pl->second.fixTopology(50);
  }

  if ((pl->second.distTo(ap) > 200) || (pl->second.distTo(bp) > 200)) {
    /**
     * something is not right, the distance from the station to its geometry
     * is excessive. fall back to straight line connection
     */
    PolyLine p = PolyLine(ap, bp);
    p.smoothenOutliers(50);
    return std::pair<bool, PolyLine>(false, p);
  }

  PolyLine p;

  if (!_cfg->ignoreGtfsDistances && distA > -1 && distA > -1 && totalTripDist > 0) {
    p = pl->second.getSegment(
        (distA - t->getShape()->getPoints().begin()->travelDist) /
            totalTripDist,
        (distB - t->getShape()->getPoints().begin()->travelDist) /
            totalTripDist);
  } else {
    p = pl->second.getSegment(ap, bp);
  }

  return std::pair<bool, PolyLine>(true, p);
}

// _____________________________________________________________________________
void Builder::averageNodePositions(Graph* g) {
  for (auto n : *g->getNodes()) {
    double x = 0;
    double y = 0;
    size_t c = 0;

    for (auto e : n->getAdjList()) {
      for (auto eg : *e->getEdgeTripGeoms()) {
        if (eg.getGeomDir() != n) {
          x += eg.getGeom().getLine().front().get<0>();
          y += eg.getGeom().getLine().front().get<1>();
        } else {
          x += eg.getGeom().getLine().back().get<0>();
          y += eg.getGeom().getLine().back().get<1>();
        }
        c++;
      }
    }

    if (c > 0) n->setPos(Point(x / c, y / c));
  }
}

// _____________________________________________________________________________
Node* Builder::addStop(gtfs::Stop* curStop, uint8_t aggrLevel, Graph* g) {
  if (aggrLevel && curStop->getParentStation() != 0) {
    return addStop(curStop->getParentStation(), aggrLevel, g);
  }

  Node* n = g->getNodeByStop(curStop, aggrLevel);

  if (n) return n;

  Point p = getProjectedPoint(curStop->getLat(), curStop->getLng(),
                              g->getProjection());

  if (aggrLevel > 1) {
    n = g->getNearestNode(p, 60);
  }

  if (n) {
    n->addStop(curStop);
  } else {
    n = g->addNode(new Node(p, curStop));
  }

  return n;
}

// _____________________________________________________________________________
bool Builder::checkTripSanity(gtfs::Trip* t) const {
  return checkShapeSanity(t->getShape());
}

// _____________________________________________________________________________
bool Builder::checkShapeSanity(gtfs::Shape* s) const {
  if (!s || s->getPoints().size() < 2) return false;
  return true;
}

// _____________________________________________________________________________
void Builder::removeEdgeArtifacts(Graph* g) {
  double MIN_SEG_LENGTH = 20;

restart:
  for (graph::Node* n : *g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (auto etg : *e->getEdgeTripGeoms()) {
        if (etg.getGeom().getLength() < MIN_SEG_LENGTH) {
          Node* from = e->getFrom();
          Node* to = e->getTo();
          if (from->getStops().size() == 0 || to->getStops().size() == 0) {
            combineNodes(from, to, g);
            // OMG really
            goto restart;
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void Builder::removeNodeArtifacts(Graph* g) {
restart:
  for (graph::Node* n : *g->getNodes()) {
    std::vector<Edge*> edges;
    edges.insert(edges.end(), n->getAdjListIn().begin(),
                 n->getAdjListIn().end());
    edges.insert(edges.end(), n->getAdjListOut().begin(),
                 n->getAdjListOut().end());
    if (edges.size() == 2 && n->getStops().size() == 0) {
      const EdgeTripGeom& etga = *edges[0]->getEdgeTripGeoms()->begin();
      const EdgeTripGeom& etgb = *edges[1]->getEdgeTripGeoms()->begin();

      if (edges[0]->getEdgeTripGeoms()->size() == 1 &&
          edges[1]->getEdgeTripGeoms()->size() == 1 &&
          etga.routeEquivalent(etgb)) {
        combineEdges(edges[0], edges[1], n, g);
        // TODO: OMG really
        goto restart;
      }
    }
  }
}

// _____________________________________________________________________________
bool Builder::combineEdges(Edge* a, Edge* b, Node* n, Graph* g) {
  assert((a->getTo() == n || a->getFrom() == n) &&
         (b->getTo() == n || b->getFrom() == n));

  if (a->getEdgeTripGeoms()->size() != 1 || b->getEdgeTripGeoms()->size() != 1)
    return false;

  EdgeTripGeom etga = *a->getEdgeTripGeoms()->begin();
  const EdgeTripGeom& etgb = *b->getEdgeTripGeoms()->begin();

  PolyLine p = etga.getGeom();

  if (etga.getGeomDir() == n) {
    if (etgb.getGeomDir() == n) {
      for (size_t i = etgb.getGeom().getLine().size(); i > 0; i--) {
        p << etgb.getGeom().getLine()[i - 1];
      }
    } else {
      for (size_t i = 0; i < etgb.getGeom().getLine().size(); i++) {
        p << etgb.getGeom().getLine()[i];
      }
    }
  } else {
    if (etgb.getGeomDir() == n) {
      for (size_t i = etgb.getGeom().getLine().size(); i > 0; i--) {
        p >> etgb.getGeom().getLine()[i - 1];
      }
    } else {
      for (size_t i = 0; i < etgb.getGeom().getLine().size(); i++) {
        p >> etgb.getGeom().getLine()[i];
      }
    }
  }

  p.smoothenOutliers(50);
  p.applyChaikinSmooth(3);
  p.simplify(3);

  Node* newFrom = a->getTo() == n ? a->getFrom() : a->getTo();
  Node* newTo = b->getFrom() == n ? b->getTo() : b->getFrom();

  if (etga.getGeomDir() == n) etga.setGeomDir(newTo);

  for (auto& to : *etga.getTripsUnordered()) {
    if (to.direction == n) to.direction = newTo;
  }

  etga.setGeom(p);

  Edge* e = g->addEdge(newFrom, newTo);
  e->addEdgeTripGeom(etga);

  g->deleteEdge(a->getFrom(), a->getTo());
  g->deleteEdge(b->getFrom(), b->getTo());
  delete (n);
  g->deleteNode(n);

  newFrom->replaceEdgeInConnections(a, e);
  newTo->replaceEdgeInConnections(b, e);

  return true;
}

// _____________________________________________________________________________
bool Builder::combineNodes(Node* a, Node* b, Graph* g) {
  assert(a->getStops().size() == 0 || b->getStops().size() == 0);

  if (a->getStops().size() != 0) {
    Node* c = a;
    a = b;
    b = c;
  }

  Edge* connecting = g->getEdge(a, b);

  std::map<const Edge*, const Edge*> oldNew;

  for (graph::Edge* e : a->getAdjListOut()) {
    if (connecting == e) continue;
    Edge* newE = g->addEdge(b, e->getTo());
    b->addEdge(newE);
    e->getTo()->removeEdge(e);

    oldNew[e] = newE;

    e->getTo()->replaceEdgeInConnections(e, newE);

    for (auto etg : *e->getEdgeTripGeoms()) {
      EdgeTripGeom etgNew(etg.getGeom(),
                          etg.getGeomDir() == a ? b : etg.getGeomDir());
      for (auto to : *etg.getTripsUnordered()) {
        if (a->isConnOccuring(to.route, e, connecting)) {
          for (const auto edge :
               b->getConnectingEdgesFor(to.route, connecting)) {
            b->connOccurs(to.route, newE, edge);
          }
        }

        for (auto trip : to.trips) {
          etgNew.addTrip(trip, to.direction == a ? b : to.direction);
        }
      }

      newE->addEdgeTripGeom(etgNew);
      newE->simplify();
    }
  }

  for (graph::Edge* e : a->getAdjListIn()) {
    if (connecting == e) continue;
    Edge* newE = g->addEdge(e->getFrom(), b);
    b->addEdge(newE);
    e->getFrom()->removeEdge(e);

    oldNew[e] = newE;

    e->getFrom()->replaceEdgeInConnections(e, newE);

    for (auto etg : *e->getEdgeTripGeoms()) {
      EdgeTripGeom etgNew(etg.getGeom(),
                          etg.getGeomDir() == a ? b : etg.getGeomDir());
      for (auto to : *etg.getTripsUnordered()) {
        if (a->isConnOccuring(to.route, e, connecting)) {
          for (const auto edge :
               b->getConnectingEdgesFor(to.route, connecting)) {
            b->connOccurs(to.route, newE, edge);
          }
        }
        for (auto trip : to.trips) {
          etgNew.addTrip(trip, to.direction == a ? b : to.direction);
        }
      }

      newE->addEdgeTripGeom(etgNew);
      newE->simplify();
    }
  }

  for (const auto& occs : a->getOccuringConnections()) {
    for (const auto& occ : occs.second) {
      if (occ.from != connecting && occ.to != connecting) {
        b->connOccurs(occs.first, oldNew[occ.from], oldNew[occ.to]);
      }
    }
  }

  g->deleteEdge(a, b);
  g->deleteNode(a);
  delete a;

  return true;
}

// _____________________________________________________________________________
PolyLine Builder::getAveragedFromSharedSeg(const ShrdSegWrap& w) const {
  const EdgeTripGeom& geomA = w.e->getEdgeTripGeoms()->front();
  const EdgeTripGeom& geomB = w.f->getEdgeTripGeoms()->front();

  PolyLine a = geomA.getGeom().getSegment(w.s.first.first.totalPos,
                                          w.s.second.first.totalPos);

  PolyLine b;

  if (w.s.first.second.totalPos > w.s.second.second.totalPos) {
    b = geomB.getGeom().getSegment(w.s.second.second.totalPos,
                                   w.s.first.second.totalPos);
    b.reverse();
  } else {
    b = geomB.getGeom().getSegment(w.s.first.second.totalPos,
                                   w.s.second.second.totalPos);
  }

  bool dominationHeuristic = false;

  if (dominationHeuristic) {
    // test if one geom is the dominant line here
    bool aDomin = lineDominatesSharedSeg(w, w.e);
    bool bDomin = lineDominatesSharedSeg(w, w.f);
    if (aDomin && !bDomin) {
      return a;
    }

    if (bDomin && !aDomin) {
      return b;
    }
  }

  std::vector<const PolyLine*> avg;
  std::vector<double> weights;
  avg.push_back(&a);
  avg.push_back(&b);
  weights.push_back(geomA.getCardinality() * geomA.getCardinality());
  weights.push_back(geomB.getCardinality() * geomB.getCardinality());

  PolyLine ret = PolyLine::average(avg, weights);
  ret.simplify(5);
  return ret;
}

// _____________________________________________________________________________
bool Builder::lineDominatesSharedSeg(const ShrdSegWrap& w, Edge* e) const {
  if (e != w.e && e != w.f) return false;

  double LOOKAHEAD = 50;
  double DELTA = .5;

  Point a;
  Point b;
  Point c;
  Point d;

  if (e == w.e) {
    const EdgeTripGeom& geom = w.e->getEdgeTripGeoms()->front();
    double lookAhead = LOOKAHEAD / geom.getGeom().getLength();
    a = geom.getGeom().getPointAt(w.s.first.first.totalPos - lookAhead).p;
    b = geom.getGeom().getPointAt(w.s.first.first.totalPos).p;
    c = geom.getGeom().getPointAt(w.s.second.first.totalPos).p;
    d = geom.getGeom().getPointAt(w.s.second.first.totalPos + lookAhead).p;
  } else if (e == w.f) {
    const EdgeTripGeom& geom = w.f->getEdgeTripGeoms()->front();
    double lookAhead = LOOKAHEAD / geom.getGeom().getLength();
    if (w.s.first.second.totalPos > w.s.second.second.totalPos) {
      a = geom.getGeom().getPointAt(w.s.first.second.totalPos + lookAhead).p;
      d = geom.getGeom().getPointAt(w.s.second.second.totalPos - lookAhead).p;
    } else {
      a = geom.getGeom().getPointAt(w.s.first.second.totalPos - lookAhead).p;
      d = geom.getGeom().getPointAt(w.s.second.second.totalPos + lookAhead).p;
    }
    b = geom.getGeom().getPointAt(w.s.first.second.totalPos).p;
    c = geom.getGeom().getPointAt(w.s.second.second.totalPos).p;
  }

  double ang = pbutil::geo::angBetween(a, b) - pbutil::geo::angBetween(c, d);
  double tang = fabs(tan(ang));

  return tang < DELTA;
}
