// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include "GraphBuilder.h"
#include "log/Log.h"
#include "./../config/TransitMapConfig.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;
using namespace gtfs;

using util::geo::Point;

// _____________________________________________________________________________
GraphBuilder::GraphBuilder(const config::Config* cfg)
: _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
}

// _____________________________________________________________________________
void GraphBuilder::consume(const Feed& f, TransitGraph* g) {
  // TODO: make this stuff configurable

  uint8_t AGGREGATE_STOPS = 2; // 1: aggregate stops already aggrgated in GTFS
                               // 2: aaggregte stops by distance

  size_t cur = 0;
  for (auto t = f.tripsBegin(); t != f.tripsEnd(); ++t) {
    cur++;
    if (t->second->getStopTimes().size() < 2) continue;

    if (t->second->getRoute()->getType() == gtfs::Route::TYPE::BUS) continue;
    if (!checkTripSanity(t->second)) continue;

    auto st = t->second->getStopTimes().begin();

    StopTime prev = *st;
    addStop(prev.getStop(), AGGREGATE_STOPS, g);
    ++st;

    for (; st != t->second->getStopTimes().end(); ++st) {
      const StopTime& cur = *st;

      Node* fromNode = g->getNodeByStop(
        prev.getStop(),
        AGGREGATE_STOPS
      );

      Node* toNode = addStop(cur.getStop(), AGGREGATE_STOPS, g);

      if (fromNode == toNode) continue;

      Edge* exE = g->getEdge(fromNode, toNode);

      if (!exE) {
        exE = g->addEdge(fromNode, toNode);
      }

      if (!exE->addTrip(t->second, toNode)) {
        std::pair<bool, geo::PolyLine> edgeGeom;
        if (AGGREGATE_STOPS) {
          Stop* frs = prev.getStop()->getParentStation() ? prev.getStop()->getParentStation() : prev.getStop();
          Stop* tos = cur.getStop()->getParentStation() ? cur.getStop()->getParentStation() : cur.getStop();
          edgeGeom = getSubPolyLine(frs, tos, t->second, prev.getShapeDistanceTravelled(), cur.getShapeDistanceTravelled(), g->getProjection());
        } else {
          edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second, prev.getShapeDistanceTravelled(), cur.getShapeDistanceTravelled(), g->getProjection());
        }

        // only take geometries that could be found using the
        // shape geometry, not some fallback
        if (edgeGeom.first) {
          exE->addTrip(t->second, edgeGeom.second, toNode,
              _cfg->lineWidth, _cfg->lineSpacing);
        }
      }
      prev = cur;
    }
  }
}

// _____________________________________________________________________________
Point GraphBuilder::getProjectedPoint(double lat, double lng, projPJ p) const {
  double x = lng;
  double y = lat;
  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;

  pj_transform(_mercProj, p, 1, 1, &x, &y, 0);

  return Point(x, y);
}

// _____________________________________________________________________________
void GraphBuilder::simplify(TransitGraph* g) {
  // try to merge both-direction edges into a single one

  for (auto n : *g->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      e->simplify();
    }
  }
}

// _____________________________________________________________________________
ShrdSegWrap GraphBuilder::getNextSharedSegment(TransitGraph* g) const {
  int i = 0;
  for (auto n : *g->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      i++;
      if (_indEdges.find(e) != _indEdges.end() || e->getEdgeTripGeoms()->size() != 1) {
        continue;
      }
      // TODO: outfactor this _______
      for (auto nt : *g->getNodes()) {
        for (auto toTest : nt->getAdjListOut()) {
          // TODO: only check edges with a SINGLE geometry atm, see also check above
          if (_indEdges.find(toTest) != _indEdges.end() || toTest->getEdgeTripGeoms()->size() != 1) {
            continue;
          }

          if (e != toTest) {
            // set dmax according to the width of the etg
            double dmax = fmax(30, e->getEdgeTripGeoms()->front().getTotalWidth() / 2 + toTest->getEdgeTripGeoms()->front().getTotalWidth() / 2 + ((e->getEdgeTripGeoms()->front().getSpacing() + toTest->getEdgeTripGeoms()->front().getSpacing()) / 2));

            geo::SharedSegments s = e->getEdgeTripGeoms()->front().getGeom().getSharedSegments(toTest->getEdgeTripGeoms()->front().getGeom(), dmax);

            if (s.segments.size() > 0) {
              _pEdges[e]++;
              if (_pEdges[e] > 20) {
                LOG(WARN) << "Too many optimiziations for " << e
                  << ", preventing further..." << std::endl;
                _indEdges.insert(e);
              }
              LOG(DEBUG) << _indEdges.size() << " / " << i << std::endl;
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
bool GraphBuilder::createTopologicalNodes(TransitGraph* g) {
  ShrdSegWrap w;
  _indEdges.clear();
  _pEdges.clear();
  bool found = false;

  while ((w = getNextSharedSegment(g)).e) {
    const EdgeTripGeom& curEdgeGeom = w.e->getEdgeTripGeoms()->front();
    const EdgeTripGeom& compEdgeGeom = w.f->getEdgeTripGeoms()->front();

    const geo::PointOnLine& eap = w.s.first.first;
    const geo::PointOnLine& ebp = w.s.second.first;
    const geo::PointOnLine& fap = w.s.first.second;
    const geo::PointOnLine& fbp = w.s.second.second;

    geo::PolyLine ea = curEdgeGeom.getGeom().getSegment(0, eap.totalPos);
    geo::PolyLine ec = curEdgeGeom.getGeom().getSegment(ebp.totalPos, 1);

    geo::PolyLine fa;
    geo::PolyLine fc;

    assert(w.f->getTo() == compEdgeGeom.getGeomDir());

    if (fap.totalPos > fbp.totalPos) {
      fa = compEdgeGeom.getGeom().getSegment(fap.totalPos, 1);
      fc = compEdgeGeom.getGeom().getSegment(0, fbp.totalPos);
    } else {
      fa = compEdgeGeom.getGeom().getSegment(0, fap.totalPos);
      fc = compEdgeGeom.getGeom().getSegment(fbp.totalPos, 1);
    }

    geo::PolyLine ab = getAveragedFromSharedSeg(w);

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

    EdgeTripGeom eaEdgeGeom(ea, a, _cfg->lineWidth, _cfg->lineSpacing);
    EdgeTripGeom abEdgeGeom(ab, b, _cfg->lineWidth, _cfg->lineSpacing);
    EdgeTripGeom ecEdgeGeom(ec, w.e->getTo(), _cfg->lineWidth, _cfg->lineSpacing);

    const Node* faDir = 0;
    const Node* fcDir = 0;

    if (fap.totalPos > fbp.totalPos) {
      faDir = w.f->getTo();
      fcDir = b;
    } else {
      faDir = a;
      fcDir = w.f->getTo();
    }

    EdgeTripGeom faEdgeGeom(fa, faDir, _cfg->lineWidth, _cfg->lineSpacing);
    EdgeTripGeom fcEdgeGeom(fc, fcDir, _cfg->lineWidth, _cfg->lineSpacing);

    for (const TripOccurance& r : curEdgeGeom.getTripsUnordered()) {
      for (auto& t : r.trips) {
        if (r.direction == w.e->getTo()) {
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
        if ((r.direction == w.f->getTo() && !(fap.totalPos > fbp.totalPos)) ||
            (r.direction != w.f->getTo() && (fap.totalPos > fbp.totalPos))) {
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

    Edge* faE = 0;
    Edge* fbE = 0;

    if (fap.totalPos > fbp.totalPos) {
      faE = g->addEdge(a, wfto);
      fbE = g->addEdge(wffrom, b);
    } else {
      faE = g->addEdge(wffrom, a);
      fbE = g->addEdge(b, wfto);
    }

    if (eaE) {
      eaE->addEdgeTripGeom(eaEdgeGeom);
      eaE->simplify();
    }
    if (abE) {
      abE->addEdgeTripGeom(abEdgeGeom);
      abE->simplify();
    }
    if (ebE) {
      ebE->addEdgeTripGeom(ecEdgeGeom);
      ebE->simplify();
    }
    if (faE) {
      faE->addEdgeTripGeom(faEdgeGeom);
      faE->simplify();
    }
    if (fbE) {
      fbE->addEdgeTripGeom(fcEdgeGeom);
      fbE->simplify();
    }

    found = true;
  }

  return found;
}

// _____________________________________________________________________________
std::pair<bool, geo::PolyLine> GraphBuilder::getSubPolyLine(Stop* a, Stop* b,
    Trip* t, double distA, double distB, projPJ proj) {
  Point ap = getProjectedPoint(a->getLat(), a->getLng(), proj);
  Point bp = getProjectedPoint(b->getLat(), b->getLng(), proj);

  if (!t->getShape()) {
    return std::pair<bool, geo::PolyLine>(false,
        geo::PolyLine(ap, bp));
  }

  double totalTripDist = t->getShape()->getPoints().rbegin()->travelDist -
    t->getShape()->getPoints().begin()->travelDist;

  auto pl = _polyLines.find(t->getShape());
  if (pl == _polyLines.end()) {
    // generate polyline for this shape
    pl = _polyLines.insert(
        std::pair<gtfs::Shape*, geo::PolyLine>(
          t->getShape(), geo::PolyLine())).first;

    for (const auto& sp : t->getShape()->getPoints()) {
      pl->second << getProjectedPoint(sp.lat, sp.lng, proj);
    }

    pl->second.simplify(20);
    pl->second.smoothenOutliers(50);
    pl->second.fixTopology(50);
  }

  if ((pl->second.distTo(ap) > 200) ||
    (pl->second.distTo(bp) > 200)) {

    /**
     * something is not right, the distance from the station to its geometry
     * is excessive. fall back to straight line connection
     */
    geo::PolyLine p = geo::PolyLine(ap, bp);
    p.smoothenOutliers(50);
    return std::pair<bool, geo::PolyLine>(false, p);
  }

  geo::PolyLine p;

  // TODO!!!!! REMOVE FALSE!!!! ONLY USE FOR FREIBURG TESTING!!!!
  // if (distA > -1 && distA > -1 && totalTripDist > 0) {
  if (distA > -1 && distA > -1 && totalTripDist > 0) {
    p = pl->second.getSegment((distA - t->getShape()->getPoints().begin()->travelDist) / totalTripDist, (distB - t->getShape()->getPoints().begin()->travelDist) / totalTripDist);
  } else {
    p = pl->second.getSegment(ap, bp);
  }

  return std::pair<bool, geo::PolyLine>(true, p);
}

// _____________________________________________________________________________
void GraphBuilder::averageNodePositions(TransitGraph* g) {
  for (auto n : *g->getNodes()) {
    double x = 0;
    double y = 0;
    size_t c = 0;

    for (auto e : n->getAdjListOut()) {
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

    for (auto e : n->getAdjListIn()) {
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

    n->setPos(Point(x/c, y/c));
  }
}

// _____________________________________________________________________________
Node* GraphBuilder::addStop(gtfs::Stop* curStop, uint8_t aggrLevel,
    TransitGraph* g) {
  if (aggrLevel && curStop->getParentStation() != 0) {
    return addStop(curStop->getParentStation(), aggrLevel, g);
  }

  Node* n = g->getNodeByStop(
    curStop,
    aggrLevel
  );

  if (n) return n;

  Point p = getProjectedPoint(curStop->getLat(), curStop->getLng(),
      g->getProjection());

  if (aggrLevel > 1) {
    n = g->getNearestNode(p, 100);
  }

  if (n) {
    n->addStop(curStop);
  } else {
    n = g->addNode(
      new Node(
        p,
        curStop
      )
    );
  }

  return n;
}

// _____________________________________________________________________________
void GraphBuilder::writeMainDirs(TransitGraph* graph) {
  for (auto n : *graph->getNodes()) {
    std::set<Edge*> eSet;
    eSet.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
    eSet.insert(n->getAdjListOut().begin(), n->getAdjListOut().end());

    for (Edge* e : eSet) {
      if (e->getEdgeTripGeoms()->size() == 0) continue;

      // atm, always take the first edge trip geometry with the most routes
      size_t lc = 0;
      const EdgeTripGeom* g = 0;
      for (auto& rg : *e->getEdgeTripGeoms()) {
         if (rg.getTripsUnordered()->size() >= lc) {
           lc = rg.getTripsUnordered()->size();
           g = &rg;
         }
      }

      if (!g || g->getGeom().getLength() == 0) continue;

      NodeFront f(e, n, g);
      geo::PolyLine pl;

      f.refEtgLengthBefExp = g->getGeom().getLength();

      if (g->getGeomDir() == n) {
        pl = g->getGeom().getOrthoLineAtDist(g->getGeom().getLength(),
            g->getTotalWidth());
      } else {
        pl = g->getGeom().getOrthoLineAtDist(0, g->getTotalWidth());
        pl.reverse();
      }

      f.setGeom(pl);

      // initial free
      freeNodeFront(&f);

      n->addMainDir(f);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::expandOverlappinFronts(TransitGraph* g) {
  // now, look at the nodes entire front geometries and expand them
  // until nothing overlaps
  double step = 1;

  while (true) {
    bool stillFree = false;
    for (auto n : *g->getNodes()) {
      std::set<NodeFront*> overlaps = nodeGetOverlappingFronts(n);
      for (auto f : overlaps) {
        stillFree = true;
        if (f->refEtg->getGeomDir() == n) {
          f->geom = f->refEtg->getGeom().getOrthoLineAtDist(
              f->refEtg->getGeom().getLength() - step,
              f->refEtg->getTotalWidth());
        } else {
          f->geom = f->refEtg->getGeom().getOrthoLineAtDist(
              step, f->refEtg->getTotalWidth());
          f->geom.reverse();
        }

        // cut the edges to fit the new front
        freeNodeFront(f);
      }
    }
    if (!stillFree) break;
  }
}

// _____________________________________________________________________________
std::set<NodeFront*> GraphBuilder::nodeGetOverlappingFronts(const Node* n)
const {
  std::set<NodeFront*> ret;
  double minLength = 1;

  for (size_t i = 0; i < n->getMainDirs().size(); ++i) {
    const NodeFront& fa = n->getMainDirs()[i];

    for (size_t j = 0; j < n->getMainDirs().size(); ++j) {
      const NodeFront& fb = n->getMainDirs()[j];

      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      if ((n->getStops().size() > 0 && fa.geom.distTo(fb.geom) < (fa.refEtg->getSpacing() + fb.refEtg->getSpacing()) / 8) ||
          (n->getStops().size() == 0 && nodeFrontsOverlap(fa, fb))) {
        if (fa.refEtg->getGeom().getLength() > minLength) {
          ret.insert(const_cast<NodeFront*>(&fa));
        }
        if (fb.refEtg->getGeom().getLength() > minLength) {
          ret.insert(const_cast<NodeFront*>(&fb));
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool GraphBuilder::nodeFrontsOverlap(const NodeFront& a,
    const NodeFront& b) const {
  size_t numShr= a.refEtg->getSharedRoutes(*b.refEtg).size();

  Point aa = a.geom.getLine().front();
  Point ab = a.geom.getLine().back();
  Point ba = b.geom.getLine().front();
  Point bb = b.geom.getLine().back();

  bool intersects = util::geo::lineIntersects(aa, ab, ba, bb);
  if (!intersects) return false;

  Point i = util::geo::intersection(aa, ab, ba, bb);

  if (numShr && a.geom.distTo(i) < (a.refEtg->getWidth() + a.refEtg->getSpacing())) return true;
  if (b.geom.distTo(a.geom) < fmin((b.refEtg->getWidth() + b.refEtg->getSpacing()) * 5, (b.refEtg->getTotalWidth() + a.refEtg->getTotalWidth()))) return true;

  return false;
}

// _____________________________________________________________________________
void GraphBuilder::freeNodeFront(NodeFront* f) {
  for (auto e : f->edges) {
    geo::PolyLine cutLine = f->geom;

    for (EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
      std::set<geo::PointOnLine, geo::PointOnLineCompare> iSects =
          cutLine.getIntersections(g.getGeom());
      if (iSects.size() > 0) {
        if (g.getGeomDir() != f->n) {
          // cut at beginning
          g.setGeom(g.getGeom().getSegment(iSects.begin()->totalPos, 1));
          assert(cutLine.distTo(g.getGeom().getLine().front()) < 0.1);
        } else {
          // cut at end
          g.setGeom(g.getGeom().getSegment(0, (--iSects.end())->totalPos));
          assert(cutLine.distTo(g.getGeom().getLine().back()) < 0.1);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::writeInitialConfig(TransitGraph* g) {
  Configuration c;
  for (graph::Node* n : *g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        Ordering order(g.getCardinality());
        for (size_t i = 0; i < g.getCardinality(); i++) order[i] = i;
        c[&g] = order;
      }
    }
  }

  g->setConfig(c);
}

// _____________________________________________________________________________
bool GraphBuilder::checkTripSanity(gtfs::Trip* t) const {
  return checkShapeSanity(t->getShape());
}

// _____________________________________________________________________________
bool GraphBuilder::checkShapeSanity(gtfs::Shape* s) const {
  if (!s|| s->getPoints().size() < 2) return false;
  return true;
}

// _____________________________________________________________________________
void GraphBuilder::removeArtifacts(TransitGraph* g) {
  double MIN_SEG_LENGTH = 20;

  restart:
  for (graph::Node* n : *g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (auto etg : *e->getEdgeTripGeoms()) {
        if (etg.getGeom().getLength() < MIN_SEG_LENGTH) {
          Node* from = e->getFrom();
          Node* to = e->getTo();
          if (from->getStops().size() == 0 || to->getStops().size() == 0) {
            g->deleteEdge(e->getFrom(), e->getTo());
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
void GraphBuilder::combineNodes(Node* a, Node* b, TransitGraph* g) {
  assert(a->getStops().size() == 0 || b->getStops().size() == 0);

  if (a->getStops().size() != 0) {
    Node* c = a;
    a = b;
    b = c;
  }

  for (graph::Edge* e : a->getAdjListOut()) {
    Edge* newE = g->addEdge(b, e->getTo());
    a->addEdge(newE);

    for (auto etg : *e->getEdgeTripGeoms()) {
      EdgeTripGeom etgNew(
          etg.getGeom(),
          etg.getGeomDir() == a ? b : etg.getGeomDir(),
          etg.getWidth(),
          etg.getSpacing());
      for (auto to : *etg.getTripsUnordered()) {
        for (auto trip : to.trips) {
          etgNew.addTrip(trip, to.direction == a ? b : to.direction);
        }
      }

      newE->addEdgeTripGeom(etgNew);
    }

    const NodeFront* nf = a->getNodeFrontFor(e);
    if (nf) {
      NodeFront newNf(newE, b, &(*newE->getEdgeTripGeoms()->begin()));
      newNf.setGeom(nf->geom);
      // NOTE: we are adding a new node front here for each edges,
      //   breaks multiple edges in one single nodefront
      b->addMainDir(newNf);
    }
  }

  for (graph::Edge* e : a->getAdjListIn()) {
    Edge* newE = g->addEdge(e->getFrom(), b);
    a->addEdge(newE);

    for (auto etg : *e->getEdgeTripGeoms()) {
      EdgeTripGeom etgNew(
          etg.getGeom(),
          etg.getGeomDir() == a ? b : etg.getGeomDir(),
          etg.getWidth(),
          etg.getSpacing());
      for (auto to : *etg.getTripsUnordered()) {
        for (auto trip : to.trips) {
          etgNew.addTrip(trip, to.direction == a ? b : to.direction);
        }
      }

      newE->addEdgeTripGeom(etgNew);
    }

    const NodeFront* nf = a->getNodeFrontFor(e);
    if (nf) {
      NodeFront newNf(newE, b, &(*newE->getEdgeTripGeoms()->begin()));
      newNf.setGeom(nf->geom);
      // NOTE: we are adding a new node front here for each edges,
      //   breaks multiple edges in one single nodefront
      b->addMainDir(newNf);
    }
  }

  g->deleteNode(a);
  delete a;
}

// _____________________________________________________________________________
geo::PolyLine GraphBuilder::getAveragedFromSharedSeg(const ShrdSegWrap& w)
const {
  const EdgeTripGeom& geomA = w.e->getEdgeTripGeoms()->front();
  const EdgeTripGeom& geomB = w.f->getEdgeTripGeoms()->front();

  geo::PolyLine a = geomA.getGeom().getSegment(w.s.first.first.totalPos,
      w.s.second.first.totalPos);
  geo::PolyLine b;

  if (w.s.first.second.totalPos > w.s.second.second.totalPos) {
    b = geomB.getGeom().getSegment(w.s.second.second.totalPos,
        w.s.first.second.totalPos);
    b.reverse();
  } else {
    b = geomB.getGeom().getSegment(w.s.first.second.totalPos,
        w.s.second.second.totalPos);
  }

  std::vector<const geo::PolyLine*> avg;
  std::vector<double> weights;
  avg.push_back(&a);
  avg.push_back(&b);
  weights.push_back(geomA.getCardinality() * geomA.getCardinality());
  weights.push_back(geomB.getCardinality() * geomB.getCardinality());

  geo::PolyLine ret = geo::PolyLine::average(avg, weights);
  ret.simplify(5);
  return ret;
}
