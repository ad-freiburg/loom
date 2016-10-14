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
GraphBuilder::GraphBuilder(TransitGraph* targetGraph, const config::Config* cfg)
: _targetGraph(targetGraph), _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
}

// _____________________________________________________________________________
void GraphBuilder::consume(const Feed& f) {
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
    addStop(prev.getStop(), AGGREGATE_STOPS);
    ++st;

    for (; st != t->second->getStopTimes().end(); ++st) {
      const StopTime& cur = *st;

      Node* fromNode = _targetGraph->getNodeByStop(
        prev.getStop(),
        AGGREGATE_STOPS
      );

      Node* toNode = addStop(cur.getStop(), AGGREGATE_STOPS);

      if (fromNode == toNode) continue;

      Edge* exE = _targetGraph->getEdge(fromNode, toNode);

      if (!exE) {
        exE = _targetGraph->addEdge(fromNode, toNode);
      }

      if (!exE->addTrip(t->second, toNode)) {
        std::pair<bool, geo::PolyLine> edgeGeom;
        if (AGGREGATE_STOPS) {
          Stop* frs = prev.getStop()->getParentStation() ? prev.getStop()->getParentStation() : prev.getStop();
          Stop* tos = cur.getStop()->getParentStation() ? cur.getStop()->getParentStation() : cur.getStop();
          edgeGeom = getSubPolyLine(frs, tos, t->second, prev.getShapeDistanceTravelled(), cur.getShapeDistanceTravelled());
        } else {
          edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second, prev.getShapeDistanceTravelled(), cur.getShapeDistanceTravelled());
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
Point GraphBuilder::getProjectedPoint(double lat, double lng) const {
  double x = lng;
  double y = lat;
  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;

  pj_transform(_mercProj, _targetGraph->getProjection(), 1, 1, &x, &y, 0);

  return Point(x, y);
}

// _____________________________________________________________________________
void GraphBuilder::simplify() {
  // try to merge both-direction edges into a single one

  for (auto n : *_targetGraph->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      e->simplify();
    }
  }
}

// _____________________________________________________________________________
ShrdSegWrap GraphBuilder::getNextSharedSegment() const {
  int i = 0;
  for (auto n : *_targetGraph->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      i++;
      if (_indEdges.find(e) != _indEdges.end() || e->getEdgeTripGeoms()->size() != 1) {
        continue;
      }
      // TODO: outfactor this _______
      for (auto nt : *_targetGraph->getNodes()) {
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
              if (_pEdges[e] > 100) {
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
bool GraphBuilder::createTopologicalNodes() {
  ShrdSegWrap w;
  _indEdges.clear();
  _pEdges.clear();
  bool found = false;

  while ((w = getNextSharedSegment()).e) {
    const EdgeTripGeom& curEdgeGeom = w.e->getEdgeTripGeoms()->front();
    const EdgeTripGeom& compEdgeGeom = w.f->getEdgeTripGeoms()->front();

    geo::PolyLine ea = curEdgeGeom.getGeom().getSegment(0, w.s.first.totalPos);
    geo::PolyLine ab = curEdgeGeom.getGeom().getSegment(w.s.first.totalPos,
        w.s.second.totalPos);
    geo::PolyLine ec = curEdgeGeom.getGeom().getSegment(w.s.second.totalPos, 1);

    geo::PointOnLine fap = compEdgeGeom.getGeom().projectOn(w.s.first.p);
    geo::PointOnLine fbp = compEdgeGeom.getGeom().projectOn(w.s.second.p);

    bool reversed = false;
    geo::PolyLine fa;
    geo::PolyLine fab;
    geo::PolyLine fc;

    assert(w.f->getTo() == compEdgeGeom.getGeomDir());

    if (fap.totalPos > fbp.totalPos) {
      reversed = true;

      fa = compEdgeGeom.getGeom().getSegment(fap.totalPos, 1);
      fab = compEdgeGeom.getGeom().getSegment(fbp, fap);
      fc = compEdgeGeom.getGeom().getSegment(0, fbp.totalPos);

      fab.reverse();
    } else {
      fa = compEdgeGeom.getGeom().getSegment(0, fap.totalPos);
      fab = compEdgeGeom.getGeom().getSegment(fap, fbp);
      fc = compEdgeGeom.getGeom().getSegment(fbp.totalPos, 1);
    }

    std::vector<const geo::PolyLine*> avg;
    avg.push_back(&ab);
    avg.push_back(&fab);
    ab = geo::PolyLine::average(avg);
    ab.simplify(5);

    //if (fab.getLength() > ab.getLength()) ab = fab;
    //

    // new nodes at the start and end of the shared segment
    Node* a = 0;
    Node* b = 0;

    double maxSnapDist = 20; //(compEdgeGeom.getTotalWidth() + curEdgeGeom.getTotalWidth()) / 2;

    if (ea.getLength() < maxSnapDist) {
      a = w.e->getFrom();
    }

    if (ec.getLength() < maxSnapDist) {
      b = w.e->getTo();
    }

    if (fa.getLength() < maxSnapDist) {
      a = !reversed ? w.f->getFrom() : w.f->getTo();
    }

    if (fc.getLength() < maxSnapDist) {
      b = !reversed ? w.f->getTo() : w.f->getFrom();
    }

    if (!a) a = new Node(w.s.first.p);
    if (!b) b = new Node(w.s.second.p);

    if (a == b) {
      continue;
    }

    EdgeTripGeom eaEdgeGeom(ea, a, _cfg->lineWidth, _cfg->lineSpacing);
    EdgeTripGeom abEdgeGeom(ab, b, _cfg->lineWidth, _cfg->lineSpacing);
    EdgeTripGeom ecEdgeGeom(ec, w.e->getTo(), _cfg->lineWidth, _cfg->lineSpacing);

    const Node* faDir = 0;
    const Node* fcDir = 0;

    if (reversed) {
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
        if ((r.direction == w.f->getTo() && !reversed) ||
            (r.direction != w.f->getTo() && reversed)) {
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
    _targetGraph->deleteEdge(w.e->getFrom(), w.e->getTo());
    _targetGraph->deleteEdge(w.f->getFrom(), w.f->getTo());

    // add new edges
    _targetGraph->addNode(a);
    _targetGraph->addNode(b);
    Edge* eaE = _targetGraph->addEdge(wefrom, a);
    Edge* abE = _targetGraph->addEdge(a, b);
    Edge* ebE = _targetGraph->addEdge(b, weto);

    Edge* faE = 0;
    Edge* fbE = 0;

    if (reversed) {
      faE =_targetGraph->addEdge(a, wfto);
      fbE =_targetGraph->addEdge(wffrom, b);
    } else {
      faE =_targetGraph->addEdge(wffrom, a);
      fbE =_targetGraph->addEdge(b, wfto);
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
    Trip* t, double distA, double distB) {
  Point ap = getProjectedPoint(a->getLat(), a->getLng());
  Point bp = getProjectedPoint(b->getLat(), b->getLng());

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
      pl->second << getProjectedPoint(sp.lat, sp.lng);
    }

    pl->second.simplify(10);
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

  if (distA > -1 && distA > -1 && totalTripDist > 0) {
    p = pl->second.getSegment((distA - t->getShape()->getPoints().begin()->travelDist) / totalTripDist, (distB - t->getShape()->getPoints().begin()->travelDist) / totalTripDist);
  } else {
    p = pl->second.getSegment(ap, bp);
  }

  return std::pair<bool, geo::PolyLine>(true, p);
}

// _____________________________________________________________________________
void GraphBuilder::averageNodePositions() {
  for (auto n : *_targetGraph->getNodes()) {
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
Node* GraphBuilder::addStop(gtfs::Stop* curStop, uint8_t aggrLevel) {
  if (aggrLevel && curStop->getParentStation() != 0) {
    return addStop(curStop->getParentStation(), aggrLevel);
  }

  Node* n = _targetGraph->getNodeByStop(
    curStop,
    aggrLevel
  );

  if (n) return n;

  Point p = getProjectedPoint(curStop->getLat(), curStop->getLng());

  if (aggrLevel > 1) {
    n = _targetGraph->getNearestNode(p, 100);
  }

  if (n) {
    n->addStop(curStop);
  } else {
    n = _targetGraph->addNode(
      new Node(
        p,
        curStop
      )
    );
  }

  return n;
}

// _____________________________________________________________________________
void GraphBuilder::writeMainDirs() {
  for (auto n : *_targetGraph->getNodes()) {
    std::set<Edge*> eSet;
    eSet.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
    eSet.insert(n->getAdjListOut().begin(), n->getAdjListOut().end());

    for (Edge* e : eSet) {
      if (e->getEdgeTripGeoms()->size() == 0) continue;

      // atm, always take the first edge trip geometry with the most routes
      size_t lc = 0;
      const EdgeTripGeom* g;
      for (auto& rg : *e->getEdgeTripGeoms()) {
         if (rg.getTripsUnordered()->size() >= lc) {
           lc = rg.getTripsUnordered()->size();
           g = &rg;
         }
      }

      if (g->getGeom().getLength() == 0) continue;

      NodeFront f(e, n, g);
      geo::PolyLine pl;

      f.refEtgLengthBefExp = g->getGeom().getLength();

      if (g->getGeomDir() == n) {
        pl = g->getGeom().getOrthoLineAtDist(g->getGeom().getLength(),
            g->getTotalWidth());
      } else {
        pl = g->getGeom().getOrthoLineAtDist(0, g->getTotalWidth());
      }

      f.setGeom(pl);

      // initial free
      freeNodeFront(&f);

      n->addMainDir(f);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::expandOverlappinFronts() {
  // now, look at the nodes entire front geometries and expand them
  // until nothing overlaps
  double step = 1;

  while (true) {
    bool stillFree = false;
    for (auto n : *_targetGraph->getNodes()) {
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

  if (numShr && a.geom.distTo(i) < a.refEtg->getWidth() + a.refEtg->getSpacing()) return true;
  if (b.geom.distTo(a.geom) < (b.refEtg->getTotalWidth() + a.refEtg->getTotalWidth()) / 2) return true;

  return a.geom.distTo(b.geom) < a.refEtg->getSpacing() + a.refEtg->getWidth();

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
void GraphBuilder::writeInitialConfig() {
  Configuration c;
  for (graph::Node* n : *_targetGraph->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        Ordering order(g.getCardinality());
        for (size_t i = 0; i < g.getCardinality(); i++) order[i] = i;
        c[&g] = order;
      }
    }
  }

  _targetGraph->setConfig(c);
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
void GraphBuilder::removeArtifacts() {
  double MIN_SEG_LENGTH = 20;

  restart:
  for (graph::Node* n : *_targetGraph->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (auto etg : *e->getEdgeTripGeoms()) {
        if (etg.getGeom().getLength() < MIN_SEG_LENGTH) {
          Node* from = e->getFrom();
          Node* to = e->getTo();
          if (from->getStops().size() == 0 || to->getStops().size() == 0) {
            _targetGraph->deleteEdge(e->getFrom(), e->getTo());
            combineNodes(from, to);
            // OMG really
            goto restart;
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::combineNodes(Node* a, Node* b) {
  assert(a->getStops().size() == 0 || b->getStops().size() == 0);

  if (a->getStops().size() != 0) {
    Node* c = a;
    a = b;
    b = c;
  }

  for (graph::Edge* e : a->getAdjListOut()) {
    Edge* newE = _targetGraph->addEdge(b, e->getTo());
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
    Edge* newE = _targetGraph->addEdge(e->getFrom(), b);
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

  _targetGraph->deleteNode(a);
  delete a;
}
