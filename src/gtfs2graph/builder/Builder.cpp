// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "gtfs2graph/builder/Builder.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "gtfs2graph/graph/EdgePL.h"
#include "gtfs2graph/graph/NodePL.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace gtfs2graph;
using namespace graph;

using util::geo::Point;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Box;
using util::geo::DBox;
using util::geo::extendBox;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using graph::Node;
using graph::Edge;

using ad::cppgtfs::gtfs::Stop;
using ad::cppgtfs::gtfs::StopTime;
using ad::cppgtfs::gtfs::Trip;
using ad::cppgtfs::gtfs::Shape;
using ad::cppgtfs::gtfs::Feed;

const static char* WGS84_PROJ =
    "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

// _____________________________________________________________________________
Builder::Builder(const config::Config* cfg) : _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
  _graphProj = pj_init_plus(cfg->projectionString.c_str());
}

// _____________________________________________________________________________
void Builder::consume(const Feed& f, BuildGraph* g) {
  DBox graphBox(getProjectedPoint(f.getMinLat(), f.getMinLon(), _graphProj),
    getProjectedPoint(f.getMaxLat(), f.getMaxLon(), _graphProj));

  NodeGrid ngrid(200, 200, graphBox);

  bool SANITY_CHECK = true;

  for (auto t = f.getTrips().begin(); t != f.getTrips().end(); ++t) {
    if (t->second->getStopTimes().size() < 2) continue;
    if (!(_cfg->useMots & (1 << t->second->getRoute()->getType()))) continue;
    if (SANITY_CHECK && !checkTripSanity(t->second)) continue;

    auto st = t->second->getStopTimes().begin();

    auto prev = *st;
    const Edge* prevEdge = 0;
    addStop(prev.getStop(), _cfg->stationAggrLevel, g, &ngrid);
    ++st;

    for (; st != t->second->getStopTimes().end(); ++st) {
      const auto& cur = *st;

      Node* fromNode = getNodeByStop(g, prev.getStop(), _cfg->stationAggrLevel);
      Node* toNode = addStop(cur.getStop(), _cfg->stationAggrLevel, g, &ngrid);

      // TODO: we should also allow this, for round-trips
      if (fromNode == toNode) continue;

      Edge* exE = g->getEdg(fromNode, toNode);

      if (!exE) {
        exE = g->addEdg(fromNode, toNode, EdgePL());
        exE->pl().setEdge(exE);
      }

      Node* directionNode = toNode;

      if (!exE->pl().addTrip(t->second, directionNode)) {
        std::pair<bool, PolyLine<double>> edgeGeom;
        if (_cfg->stationAggrLevel) {
          const Stop* frs = prev.getStop()->getParentStation()
                                ? prev.getStop()->getParentStation()
                                : prev.getStop();
          const Stop* tos = cur.getStop()->getParentStation()
                                ? cur.getStop()->getParentStation()
                                : cur.getStop();
          edgeGeom = getSubPolyLine(frs, tos, t->second,
                                    prev.getShapeDistanceTravelled(),
                                    cur.getShapeDistanceTravelled());
        } else {
          edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second,
                                    prev.getShapeDistanceTravelled(),
                                    cur.getShapeDistanceTravelled());
        }

        // only take geometries that could be found using the
        // shape geometry, not some fallback
        if (!SANITY_CHECK || edgeGeom.first) {
          if (prevEdge) {
            fromNode->pl().connOccurs(t->second->getRoute(), prevEdge, exE);
          }

          exE->pl().addTrip(t->second, edgeGeom.second, directionNode);
        }
      }
      prev = cur;
      prevEdge = exE;
    }
  }
}

// _____________________________________________________________________________
DPoint Builder::getProjectedPoint(double lat, double lng, projPJ p) const {
  double x = lng * DEG_TO_RAD, y = lat * DEG_TO_RAD;
  pj_transform(_mercProj, p, 1, 1, &x, &y, 0);
  return DPoint(x, y);
}

// _____________________________________________________________________________
void Builder::simplify(BuildGraph* g) {
  // delete edges without a reference ETG
  std::vector<Edge*> toDel;
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (!e->pl().getRefETG()) toDel.push_back(e);
    }
  }

  for (auto e : toDel) g->delEdg(e->getFrom(), e->getTo());

  // try to merge both-direction edges into a single one
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      e->pl().simplify();
    }
  }
}

// _____________________________________________________________________________
ShrdSegWrap Builder::getNextSharedSegment(BuildGraph* g, bool final, EdgeGrid* grid) const {
  int i = 0;
  double dmin = _cfg->maxAggrDistance;
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      i++;
      if (_indEdges.find(e) != _indEdges.end() ||
          e->pl().getEdgeTripGeoms()->size() != 1) {
        continue;
      }

      std::set<Edge*> neighbors;
      grid->getNeighbors(e, fmax(5, dmin * 10), &neighbors);

      for (auto toTest : neighbors) {
        // TODO: only check edges with a SINGLE geometry atm, see also check
        // above
        if (_indEdgesPairs.find(std::pair<const Edge*, const Edge*>(
                e, toTest)) != _indEdgesPairs.end() ||
            _indEdges.find(toTest) != _indEdges.end() ||
            toTest->pl().getEdgeTripGeoms()->size() != 1) {
          continue;
        }

        if (e != toTest) {
          double dmax = 5;

          if (final) {
            dmax =
                fmax(dmin, e->pl().getRefETG()->getCardinality() * (dmin / 2) +
                               toTest->pl().getRefETG()->getCardinality() *
                                   (dmin / 2));
          }

          SharedSegments<double> s = e->pl().getRefETG()->getGeom().getSharedSegments(
              toTest->pl().getRefETG()->getGeom(), dmax);

          if (s.segments.size() > 0) {
            _pEdges[std::pair<const Edge*, const Edge*>(e, toTest)]++;
            _pEdges[std::pair<const Edge*, const Edge*>(toTest, e)]++;

            if (_pEdges[std::pair<const Edge*, const Edge*>(e, toTest)] > 20) {
              _indEdgesPairs.insert(
                  std::pair<const Edge*, const Edge*>(e, toTest));
              _indEdgesPairs.insert(
                  std::pair<const Edge*, const Edge*>(toTest, e));
            }

            return ShrdSegWrap(e, toTest, s.segments.front());
          }
        }
      }

      // nothing overlapping found for this edge anymore, mark as processed
      _indEdges.insert(e);
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

  double lookDist = 100, minInAngle = 0.2;

  const EdgeTripGeom& e = *ea->pl().getRefETG();
  const EdgeTripGeom& f = *eb->pl().getRefETG();

  double angleE = 0, angleF = 0, angleEPartner = 0, angleFPartner = 0;

  if (e.getGeomDir() == a) {
    // look at end of geom
    DPoint a = e.getGeom().back();
    DPoint b = e.getGeom().getPointAtDist(e.getGeom().getLength() - lookDist).p;
    angleE = util::geo::angBetween(a, b);
  } else {
    // look at front of geom
    DPoint a = e.getGeom().front();
    DPoint b = e.getGeom().getPointAtDist(lookDist).p;
    angleE = util::geo::angBetween(a, b);
  }

  if (f.getGeomDir() == a) {
    // look at end of geom
    DPoint a = f.getGeom().back();
    DPoint b = f.getGeom().getPointAtDist(f.getGeom().getLength() - lookDist).p;
    angleF = util::geo::angBetween(a, b);
  } else {
    // look at front of geom
    DPoint a = f.getGeom().front();
    DPoint b = f.getGeom().getPointAtDist(lookDist).p;
    angleF = util::geo::angBetween(a, b);
  }

  double inDiff = fabs(angleE - angleF);

  if (inDiff < minInAngle) return false;

  bool eFound = false, fFound = false;

  // now lets get some edge that is not e or f and has all of e's routes

  for (auto edge : a->getAdjList()) {
    if (edge == ea || edge == eb) continue;
    if (edge->pl().getEdgeTripGeoms()->size() == 0) continue;
    const EdgeTripGeom& etg = *edge->pl().getRefETG();

    if (etg.routeEquivalent(e)) {
      if (etg.getGeomDir() != a) {
        // look at end of geom
        DPoint a = etg.getGeom().back();
        DPoint b = etg.getGeom()
                      .getPointAtDist(etg.getGeom().getLength() - lookDist)
                      .p;
        angleEPartner = util::geo::angBetween(a, b);
      } else {
        // look at front of geom
        DPoint a = etg.getGeom().front();
        DPoint b = etg.getGeom().getPointAtDist(lookDist).p;
        angleEPartner = util::geo::angBetween(a, b);
      }
      eFound = true;
      break;
    }
  }

  for (auto edge : a->getAdjList()) {
    if (edge == ea || edge == eb) continue;
    if (edge->pl().getEdgeTripGeoms()->size() == 0) continue;
    const EdgeTripGeom& etg = *edge->pl().getRefETG();

    if (etg.routeEquivalent(f)) {
      if (etg.getGeomDir() != a) {
        // look at end of geom
        DPoint a = etg.getGeom().back();
        DPoint b = etg.getGeom()
                      .getPointAtDist(etg.getGeom().getLength() - lookDist)
                      .p;
        angleFPartner = util::geo::angBetween(a, b);
      } else {
        // look at front of geom
        DPoint a = etg.getGeom().front();
        DPoint b = etg.getGeom().getPointAtDist(lookDist).p;
        angleFPartner = util::geo::angBetween(a, b);
      }
      fFound = true;
      break;
    }
  }

  if (!eFound || !fFound) return false;

  double outDiff = fabs(angleEPartner - angleFPartner);

  if (outDiff < minInAngle) return false;

  if (fabs(angleE - angleEPartner) < 0.2 && fabs(angleF - angleFPartner) < 0.2) {
    return true;
  }

  return false;
}

// _____________________________________________________________________________
bool Builder::createTopologicalNodes(BuildGraph* g, bool final) {
  ShrdSegWrap w;
  _indEdges.clear();
  _pEdges.clear();
  bool found = false;

  EdgeGrid grid = getGeoIndex(g);

  while ((w = getNextSharedSegment(g, final, &grid)).e) {
    const EdgeTripGeom& curEdgeGeom = *w.e->pl().getRefETG();
    const EdgeTripGeom& cmpEdgeGeom = *w.f->pl().getRefETG();

    const LinePoint<double>& eap = w.s.first.first;
    const LinePoint<double>& ebp = w.s.second.first;
    const LinePoint<double>& fap = w.s.first.second;
    const LinePoint<double>& fbp = w.s.second.second;

    PolyLine<double> ea = curEdgeGeom.getGeom().getSegment(0, eap.totalPos);
    PolyLine<double> ec = curEdgeGeom.getGeom().getSegment(ebp.totalPos, 1);

    PolyLine<double> fa, fc;

    assert(w.f->getTo() == cmpEdgeGeom.getGeomDir());

    if (fap.totalPos > fbp.totalPos) {
      fa = cmpEdgeGeom.getGeom().getSegment(fap.totalPos, 1);
      fc = cmpEdgeGeom.getGeom().getSegment(0, fbp.totalPos);
    } else {
      fa = cmpEdgeGeom.getGeom().getSegment(0, fap.totalPos);
      fc = cmpEdgeGeom.getGeom().getSegment(fbp.totalPos, 1);
    }

    PolyLine<double> ab = getAveragedFromSharedSeg(w);

    // new nodes at the start and end of the shared segment
    Node *a = 0, *b = 0;

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

    if (!a) {
      a = g->addNd(NodePL(w.s.first.first.p));
      a->pl().setNode(a);
    }
    if (!b) {
      b = g->addNd(NodePL(w.s.second.first.p));
      b->pl().setNode(b);
    }

    if (a == b) {
      continue;
    }

    // TODO: what about the nodes we added above??
    if (lineCrossesAtNode(a, w.e, w.f)) continue;
    if (lineCrossesAtNode(b, w.e, w.f)) continue;

    EdgeTripGeom eaEdgeGeom(ea, a);
    EdgeTripGeom abEdgeGeom(ab, b);
    EdgeTripGeom ecEdgeGeom(ec, w.e->getTo());

    const Node *faDir = 0, *fcDir = 0;

    if (fap.totalPos > fbp.totalPos) {
      faDir = w.f->getTo();
      fcDir = b;
    } else {
      faDir = a;
      fcDir = w.f->getTo();
    }

    EdgeTripGeom faEdgeGeom(fa, faDir);
    EdgeTripGeom fcEdgeGeom(fc, fcDir);

    Node* wefrom = w.e->getFrom();
    Node* weto = w.e->getTo();
    Node* wffrom = w.f->getFrom();
    Node* wfto = w.f->getTo();

    for (const TripOccurance& r : curEdgeGeom.getTripsUnordered()) {
      for (auto& t : r.trips) {
        if (!r.direction) {
          eaEdgeGeom.addTrip(t, 0);
          abEdgeGeom.addTrip(t, 0);
          ecEdgeGeom.addTrip(t, 0);
        } else if (r.direction == weto) {
          eaEdgeGeom.addTrip(t, a);
          abEdgeGeom.addTrip(t, b);
          ecEdgeGeom.addTrip(t, weto);
        } else {
          eaEdgeGeom.addTrip(t, wefrom);
          abEdgeGeom.addTrip(t, a);
          ecEdgeGeom.addTrip(t, b);
        }
      }
    }

    for (const TripOccurance r : cmpEdgeGeom.getTripsUnordered()) {
      for (auto& t : r.trips) {
        if (!r.direction) {
          faEdgeGeom.addTrip(t, 0);
          abEdgeGeom.addTrip(t, 0);
          fcEdgeGeom.addTrip(t, 0);
        } else if ((r.direction == wfto)) {
          if (fap.totalPos > fbp.totalPos) {
            faEdgeGeom.addTrip(t, wfto);
            abEdgeGeom.addTrip(t, a);
            fcEdgeGeom.addTrip(t, b);
          } else {
            faEdgeGeom.addTrip(t, a);
            abEdgeGeom.addTrip(t, b);
            fcEdgeGeom.addTrip(t, wfto);
          }
        } else {
          if (fap.totalPos > fbp.totalPos) {
            faEdgeGeom.addTrip(t, a);
            abEdgeGeom.addTrip(t, b);
            fcEdgeGeom.addTrip(t, b);
          } else {
            faEdgeGeom.addTrip(t, wffrom);
            abEdgeGeom.addTrip(t, a);
            fcEdgeGeom.addTrip(t, wfto);
          }
        }
      }
    }

    for (auto to : curEdgeGeom.getTripsUnordered()) {
      assert(abEdgeGeom.containsRoute(to.route));
    }

    for (auto to : cmpEdgeGeom.getTripsUnordered()) {
      assert(abEdgeGeom.containsRoute(to.route));
    }

    // delete old edges
    grid.remove(g->getEdg(w.e->getFrom(), w.e->getTo()));
    grid.remove(g->getEdg(w.f->getFrom(), w.f->getTo()));
    g->delEdg(w.e->getFrom(), w.e->getTo());
    g->delEdg(w.f->getFrom(), w.f->getTo());

    Edge *eaE = 0, *abE = 0, *ebE = 0, *faE = 0, *fbE = 0;

    // add new edges
    assert(a != b);
    abE = g->addEdg(a, b, EdgePL());
    if (abE) abE->pl().setEdge(abE);

    if (a != wefrom) {
      eaE = g->addEdg(wefrom, a, EdgePL());
      eaE->pl().setEdge(eaE);
      wefrom->pl().replaceEdgeInConnections(w.e, eaE);
    } else {
      a->pl().replaceEdgeInConnections(w.e, abE);
    }

    if (b != weto) {
      ebE = g->addEdg(b, weto, EdgePL());
      ebE->pl().setEdge(ebE);
      weto->pl().replaceEdgeInConnections(w.e, ebE);
    } else {
      assert(b == weto);
      b->pl().replaceEdgeInConnections(w.e, abE);
    }


    if (fap.totalPos > fbp.totalPos) {
      if (a != wfto) {
        faE = g->addEdg(a, wfto, EdgePL());
        faE->pl().setEdge(faE);
        wfto->pl().replaceEdgeInConnections(w.f, faE);
      } else {
        assert(a == wfto);
        a->pl().replaceEdgeInConnections(w.f, abE);
      }

      if (b != wffrom) {
        fbE = g->addEdg(wffrom, b, EdgePL());
        if (fbE) fbE->pl().setEdge(fbE);
        wffrom->pl().replaceEdgeInConnections(w.f, fbE);
      } else {
        b->pl().replaceEdgeInConnections(w.f, abE);
      }
    } else {

      if (a != wffrom) {
        faE = g->addEdg(wffrom, a, EdgePL());
        faE->pl().setEdge(faE);
        wffrom->pl().replaceEdgeInConnections(w.f, faE);
      } else {
        assert(a == wffrom);
        a->pl().replaceEdgeInConnections(w.f, abE);
      }

      if (b != wfto) {
        fbE = g->addEdg(b, wfto, EdgePL());
        fbE->pl().setEdge(fbE);
        wfto->pl().replaceEdgeInConnections(w.f, fbE);
      } else {
        b->pl().replaceEdgeInConnections(w.f, abE);
      }
    }

    if (abE) {
      abE->pl().addEdgeTripGeom(abEdgeGeom);
      abE->pl().simplify();
      grid.add(abE->pl().getRefETG()->getGeom().getLine(), abE);
    } else {
      // we use abE below without checking!
      assert(false);
    }

    if (eaE) {
      eaE->pl().addEdgeTripGeom(eaEdgeGeom);
      eaE->pl().simplify();
      grid.add(eaE->pl().getRefETG()->getGeom().getLine(), eaE);
      a->pl().sewConnectionsTogether(eaE, abE);
    }
    if (ebE) {
      ebE->pl().addEdgeTripGeom(ecEdgeGeom);
      ebE->pl().simplify();
      grid.add(ebE->pl().getRefETG()->getGeom().getLine(), ebE);
      b->pl().sewConnectionsTogether(abE, ebE);
    }

    if (faE) {
      faE->pl().addEdgeTripGeom(faEdgeGeom);
      faE->pl().simplify();
      grid.add(faE->pl().getRefETG()->getGeom().getLine(), faE);
      a->pl().sewConnectionsTogether(faE, abE);
    }
    if (fbE) {
      fbE->pl().addEdgeTripGeom(fcEdgeGeom);
      fbE->pl().simplify();
      grid.add(fbE->pl().getRefETG()->getGeom().getLine(), fbE);
      b->pl().sewConnectionsTogether(abE, fbE);
    }

    found = true;
  }

  return found;
}

// _____________________________________________________________________________
std::pair<bool, PolyLine<double>> Builder::getSubPolyLine(const Stop* a, const Stop* b,
                                                  Trip* t, double distA,
                                                  double distB) {
  DPoint ap = getProjectedPoint(a->getLat(), a->getLng(), _graphProj);
  DPoint bp = getProjectedPoint(b->getLat(), b->getLng(), _graphProj);

  if (!t->getShape()) {
    return std::pair<bool, PolyLine<double>>(false, PolyLine<double>(ap, bp));
  }

  double totalTripDist = t->getShape()->getPoints().rbegin()->travelDist -
                         t->getShape()->getPoints().begin()->travelDist;

  auto pl = _polyLines.find(t->getShape());
  if (pl == _polyLines.end()) {
    // generate polyline for this shape
    pl = _polyLines
             .insert(std::pair<Shape*, PolyLine<double>>(t->getShape(), PolyLine<double>()))
             .first;

    for (const auto& sp : t->getShape()->getPoints()) {
      pl->second << getProjectedPoint(sp.lat, sp.lng, _graphProj);
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
    PolyLine<double> p = PolyLine<double>(ap, bp);
    p.smoothenOutliers(50);
    return std::pair<bool, PolyLine<double>>(false, p);
  }

  PolyLine<double> p;

  if (!_cfg->ignoreGtfsDistances && distA > -1 && distA > -1 &&
      totalTripDist > 0) {
    p = pl->second.getSegment(
        (distA - t->getShape()->getPoints().begin()->travelDist) /
            totalTripDist,
        (distB - t->getShape()->getPoints().begin()->travelDist) /
            totalTripDist);
  } else {
    p = pl->second.getSegment(ap, bp);
  }

  return std::pair<bool, PolyLine<double>>(true, p);
}

// _____________________________________________________________________________
void Builder::averageNodePositions(BuildGraph* g) {
  for (auto n : *g->getNds()) {
    double x = 0, y = 0;
    size_t c = 0;

    for (auto e : n->getAdjList()) {
      for (auto eg : *e->pl().getEdgeTripGeoms()) {
        if (eg.getGeomDir() != n) {
          x += eg.getGeom().front().getX();
          y += eg.getGeom().front().getY();
        } else {
          x += eg.getGeom().back().getX();
          y += eg.getGeom().back().getY();
        }
        c++;
      }
    }

    if (c > 0) n->pl().setPos(DPoint(x / static_cast<double>(c), y / static_cast<double>(c)));
  }
}

// _____________________________________________________________________________
Node* Builder::addStop(const Stop* curStop, uint8_t aggrLevel, BuildGraph* g, NodeGrid* grid) {
  if (aggrLevel && curStop->getParentStation() != 0) {
    Node* n = addStop(curStop->getParentStation(), aggrLevel, g, grid);
    _stopNodes[curStop] = n;
    return n;
  }

  Node* n = getNodeByStop(g, curStop, aggrLevel);
  if (n) return n;

  DPoint p = getProjectedPoint(curStop->getLat(), curStop->getLng(), _graphProj);

  if (aggrLevel > 1) {
    n = getNearestStop(g, p, _cfg->stationAggrDistance, grid);
  }

  if (n) {
    n->pl().addStop(curStop);
    _stopNodes[curStop] = n;
  } else {
    n = g->addNd(NodePL(p, curStop));
    n->pl().setNode(n);
    grid->add(n->pl().getPos(), n);
    _stopNodes[curStop] = n;
  }

  return n;
}

// _____________________________________________________________________________
bool Builder::checkTripSanity(Trip* t) const {
  return checkShapeSanity(t->getShape());
}

// _____________________________________________________________________________
bool Builder::checkShapeSanity(Shape* s) const {
  if (!s || s->getPoints().size() < 2) return false;
  return true;
}

// _____________________________________________________________________________
void Builder::removeEdgeArtifacts(BuildGraph* g) {
  double MIN_SEG_LENGTH = 20;

restart:
  for (Node* n : *g->getNds()) {
    for (Edge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto etg : *e->pl().getEdgeTripGeoms()) {
        if (etg.getGeom().getLength() < MIN_SEG_LENGTH) {
          Node* from = e->getFrom();
          Node* to = e->getTo();
          if (from->pl().getStops().size() == 0 ||
              to->pl().getStops().size() == 0) {
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
void Builder::removeNodeArtifacts(BuildGraph* g) {
restart:
  for (Node* n : *g->getNds()) {
    std::vector<Edge*> edges;
    edges.insert(edges.end(), n->getAdjList().begin(),
                 n->getAdjList().end());
    if (edges.size() == 2 && n->pl().getStops().size() == 0) {
      const EdgeTripGeom& etga = *edges[0]->pl().getEdgeTripGeoms()->begin();
      const EdgeTripGeom& etgb = *edges[1]->pl().getEdgeTripGeoms()->begin();

      if (edges[0]->pl().getEdgeTripGeoms()->size() == 1 &&
          edges[1]->pl().getEdgeTripGeoms()->size() == 1 &&
          etga.routeEquivalent(etgb)) {
        combineEdges(edges[0], edges[1], n, g);
        // TODO: OMG really
        goto restart;
      }
    }
  }
}

// _____________________________________________________________________________
bool Builder::combineEdges(Edge* a, Edge* b, Node* n, BuildGraph* g) {
  assert((a->getTo() == n || a->getFrom() == n) &&
         (b->getTo() == n || b->getFrom() == n));

  if (a->pl().getEdgeTripGeoms()->size() != 1 ||
      b->pl().getEdgeTripGeoms()->size() != 1)
    return false;

  EdgeTripGeom etga = *a->pl().getEdgeTripGeoms()->begin();
  const EdgeTripGeom& etgb = *b->pl().getEdgeTripGeoms()->begin();

  PolyLine<double> p = etga.getGeom();

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

  assert(newFrom != newTo);

  Edge* e = g->addEdg(newFrom, newTo, EdgePL());
  e->pl().setEdge(e);
  e->pl().addEdgeTripGeom(etga);

  g->delEdg(a->getFrom(), a->getTo());
  g->delEdg(b->getFrom(), b->getTo());
  g->delNd(n);

  newFrom->pl().replaceEdgeInConnections(a, e);
  newTo->pl().replaceEdgeInConnections(b, e);

  return true;
}

// _____________________________________________________________________________
bool Builder::combineNodes(Node* a, Node* b, BuildGraph* g) {
  assert(a->pl().getStops().size() == 0 || b->pl().getStops().size() == 0);

  if (a->pl().getStops().size() != 0) {
    Node* c = a;
    a = b;
    b = c;
  }

  Edge* connecting = g->getEdg(a, b);

  std::map<const Edge*, const Edge*> oldNew;

  for (Edge* e : a->getAdjList()) {
    if (e->getFrom() != a) continue;
    if (connecting == e) continue;
    assert(b != e->getTo());
    Edge* newE = g->addEdg(b, e->getTo(), EdgePL());
    newE->pl().setEdge(newE);
    b->addEdge(newE);
    e->getTo()->removeEdge(e);

    oldNew[e] = newE;

    e->getTo()->pl().replaceEdgeInConnections(e, newE);

    for (auto etg : *e->pl().getEdgeTripGeoms()) {
      EdgeTripGeom etgNew(etg.getGeom(),
                          etg.getGeomDir() == a ? b : etg.getGeomDir());
      for (auto to : *etg.getTripsUnordered()) {
        if (a->pl().isConnOccuring(to.route, e, connecting)) {
          for (const auto edge :
               b->pl().getConnectingEdgesFor(to.route, connecting)) {
            b->pl().connOccurs(to.route, newE, edge);
          }
        }

        for (auto trip : to.trips) {
          etgNew.addTrip(trip, to.direction == a ? b : to.direction);
        }
      }

      newE->pl().addEdgeTripGeom(etgNew);
      newE->pl().simplify();
    }
  }

  for (Edge* e : a->getAdjListIn()) {
    if (e->getTo() != a) continue;
    if (connecting == e) continue;
    assert(b != e->getFrom());
    Edge* newE = g->addEdg(e->getFrom(), b, EdgePL());
    newE->pl().setEdge(newE);
    b->addEdge(newE);
    e->getFrom()->removeEdge(e);

    oldNew[e] = newE;

    e->getFrom()->pl().replaceEdgeInConnections(e, newE);

    for (auto etg : *e->pl().getEdgeTripGeoms()) {
      EdgeTripGeom etgNew(etg.getGeom(),
                          etg.getGeomDir() == a ? b : etg.getGeomDir());
      for (auto to : *etg.getTripsUnordered()) {
        if (a->pl().isConnOccuring(to.route, e, connecting)) {
          for (const auto edge :
               b->pl().getConnectingEdgesFor(to.route, connecting)) {
            b->pl().connOccurs(to.route, newE, edge);
          }
        }
        for (auto trip : to.trips) {
          etgNew.addTrip(trip, to.direction == a ? b : to.direction);
        }
      }

      newE->pl().addEdgeTripGeom(etgNew);
      newE->pl().simplify();
    }
  }

  for (const auto& occs : a->pl().getOccuringConnections()) {
    for (const auto& occ : occs.second) {
      if (occ.from != connecting && occ.to != connecting) {
        b->pl().connOccurs(occs.first, oldNew[occ.from], oldNew[occ.to]);
      }
    }
  }

  g->delEdg(a, b);
  g->delNd(a);

  return true;
}

// _____________________________________________________________________________
PolyLine<double> Builder::getAveragedFromSharedSeg(const ShrdSegWrap& w) const {
  const EdgeTripGeom& geomA = *w.e->pl().getRefETG();
  const EdgeTripGeom& geomB = *w.f->pl().getRefETG();

  PolyLine<double> a = geomA.getGeom().getSegment(w.s.first.first.totalPos,
                                          w.s.second.first.totalPos);

  PolyLine<double> b;

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

  std::vector<const PolyLine<double>*> avg{&a, &b};
  std::vector<double> weights{
      1.0 * geomA.getCardinality() * geomA.getCardinality(),
      1.0 * geomB.getCardinality() * geomB.getCardinality()};

  PolyLine<double> ret = PolyLine<double>::average(avg, weights);
  ret.simplify(5);
  return ret;
}

// _____________________________________________________________________________
bool Builder::lineDominatesSharedSeg(const ShrdSegWrap& w, Edge* e) const {
  if (e != w.e && e != w.f) return false;

  double LOOKAHEAD = 50, DELTA = .5;

  DPoint a, b, c, d;

  if (e == w.e) {
    const EdgeTripGeom& geom = *w.e->pl().getRefETG();
    double lookAhead = LOOKAHEAD / geom.getGeom().getLength();
    a = geom.getGeom().getPointAt(w.s.first.first.totalPos - lookAhead).p;
    b = geom.getGeom().getPointAt(w.s.first.first.totalPos).p;
    c = geom.getGeom().getPointAt(w.s.second.first.totalPos).p;
    d = geom.getGeom().getPointAt(w.s.second.first.totalPos + lookAhead).p;
  } else if (e == w.f) {
    const EdgeTripGeom& geom = *w.f->pl().getRefETG();
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

  double ang = util::geo::angBetween(a, b) - util::geo::angBetween(c, d);
  double tang = fabs(tan(ang));

  return tang < DELTA;
}

// _____________________________________________________________________________
Node* Builder::getNodeByStop(const BuildGraph* g, const gtfs::Stop* s,
                             bool getParent) const {
  if (getParent && s->getParentStation())
    return getNodeByStop(g, s->getParentStation());

  return getNodeByStop(g, s);
}

// _____________________________________________________________________________
Node* Builder::getNodeByStop(const BuildGraph* g, const gtfs::Stop* s) const {
  if (_stopNodes.find(s) != _stopNodes.end()) return _stopNodes.find(s)->second;

  for (const auto n : g->getNds()) {
    if (n->pl().getStops().find(const_cast<gtfs::Stop*>(s)) !=
        n->pl().getStops().end()) {
      return n;
    }
  }
  return 0;
}

// _____________________________________________________________________________
Node* Builder::getNearestStop(const BuildGraph* g, const DPoint& p,
                              double maxD, const NodeGrid* grid) const {
  double curD = DBL_MAX;

  Node* curN = 0;
  std::set<Node*> neighbors;
  grid->get(p, maxD, &neighbors);

  for (auto n : neighbors) {
    double d = util::geo::dist(n->pl().getPos(), p);
    if (d < maxD && d < curD) {
      curN = n;
      curD = d;
    }
  }

  return curN;
}

// _____________________________________________________________________________
EdgeGrid Builder::getGeoIndex(const BuildGraph* g) const {
  EdgeGrid grid(120, 120, getGraphBoundingBox(g));

  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (!e->pl().getRefETG()) continue;
      grid.add(e->pl().getRefETG()->getGeom().getLine(), e);
    }
  }

  return grid;
}

// _____________________________________________________________________________
DBox Builder::getGraphBoundingBox(const BuildGraph* g) const {
  DBox b;

  for (auto n : g->getNds()) {
    b = extendBox(n->pl().getPos(), b);
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (!e->pl().getRefETG()) {
        continue;
      }
      b = extendBox(e->pl().getRefETG()->getGeom().getLine(), b);
    }
  }

  return b;
}
