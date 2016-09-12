// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <proj_api.h>
#include "GraphBuilder.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;
using namespace gtfs;

// _____________________________________________________________________________ 
GraphBuilder::GraphBuilder(TransitGraph* targetGraph)
: _targetGraph(targetGraph) {
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
    if (!(t->second->getShape())) continue;
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
        geo::PolyLine edgeGeom;
        if (AGGREGATE_STOPS) {
          Stop* frs = prev.getStop()->getParentStation() ? prev.getStop()->getParentStation() : prev.getStop();
          Stop* tos = cur.getStop()->getParentStation() ? cur.getStop()->getParentStation() : cur.getStop();
          edgeGeom = getSubPolyLine(frs, tos, t->second);
        } else {
          edgeGeom = getSubPolyLine(prev.getStop(), cur.getStop(), t->second);
        }
        exE->addTrip(t->second, edgeGeom, toNode);
      }
      prev = cur;
    }
  }
}

// _____________________________________________________________________________ 
util::geo::Point GraphBuilder::getProjectedPoint(double lat, double lng) const {
  double x = lng;
  double y = lat;
  x *= DEG_TO_RAD;
  y *= DEG_TO_RAD;

  pj_transform(_mercProj, _targetGraph->getProjection(), 1, 1, &x, &y, 0);

  return util::geo::Point(x, y);
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
      if (_indEdges.find(e) != _indEdges.end() || e->getEdgeTripGeoms()->size() != 1) continue;
      // TODO: outfactor this _______
      for (auto nt : *_targetGraph->getNodes()) {
        for (auto toTest : nt->getAdjListOut()) {
          // TODO: only check edges with a SINGLE geometry atm, see also check above
          if (_indEdges.find(toTest) != _indEdges.end() || toTest->getEdgeTripGeoms()->size() != 1) continue;
          if (e != toTest) {
            // set dmax according to the width of the etg
            double dmax = fmax(30, e->getEdgeTripGeoms()->front().getTotalWidth() / 2 + toTest->getEdgeTripGeoms()->front().getTotalWidth() / 2);
            geo::SharedSegments s = e->getEdgeTripGeoms()->front().getGeom().getSharedSegments(toTest->getEdgeTripGeoms()->front().getGeom(), dmax);

            if (s.segments.size() > 0 && util::geo::dist(s.segments[0].first.p, s.segments[0].second.p) > 50) {
              _pEdges[e]++;
              if (_pEdges[e] > 100) {
                std::cout << "Too many optimiziations for " << e << ", preventing further..." << std::endl;
                _indEdges.insert(e);
              }
              std::cout << _indEdges.size() << " / " << i << std::endl;
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
void GraphBuilder::createTopologicalNodes() {
  ShrdSegWrap w;
  while ((w = getNextSharedSegment()).e) {
    const EdgeTripGeom& curEdgeGeom = w.e->getEdgeTripGeoms()->front();
    const EdgeTripGeom& compEdgeGeom = w.f->getEdgeTripGeoms()->front();

    geo::PolyLine ea = curEdgeGeom.getGeom().getSegment(0, w.s.first.totalPos);
    geo::PolyLine ab = curEdgeGeom.getGeom().getSegment(w.s.first.totalPos, w.s.second.totalPos);
    geo::PolyLine ec = curEdgeGeom.getGeom().getSegment(w.s.second.totalPos, 1);

    geo::PointOnLine fap = compEdgeGeom.getGeom().projectOn(w.s.first.p);
    geo::PointOnLine fbp = compEdgeGeom.getGeom().projectOn(w.s.second.p);

    bool reversed = false;
    geo::PolyLine fa;
    geo::PolyLine fab;
    geo::PolyLine fc;
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
      a = !reversed ? w.f->getFrom() : w.f->getTo();
    }

    if (fc.getLength() < maxSnapDist) {
      b = !reversed ? w.f->getTo() : w.f->getFrom();
    }

    if (!a) a = new Node(w.s.first.p);
    if (!b) b = new Node(w.s.second.p);

    EdgeTripGeom eaEdgeGeom(ea, a);
    EdgeTripGeom abEdgeGeom(ab, b);
    EdgeTripGeom ecEdgeGeom(ec, w.e->getTo());

    const Node* faDir = 0;
    const Node* fcDir = 0;

    if (reversed) {
      faDir = w.f->getTo();
      fcDir = b;
    } else {
      faDir = a;
      fcDir = w.f->getTo();
    }

    EdgeTripGeom faEdgeGeom(fa, faDir);
    EdgeTripGeom fcEdgeGeom(fc, fcDir);

    for (size_t i : curEdgeGeom.getTripOrdering()) {
      const TripOccurance& r = curEdgeGeom.getTripsUnordered()[i];
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

    for (size_t i : compEdgeGeom.getTripOrdering()) {
      const TripOccurance& r = compEdgeGeom.getTripsUnordered()[i];
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

    Node* wefrom = w.e->getFrom();
    Node* weto = w.e->getTo();
    Node* wffrom = w.f->getFrom();
    Node* wfto = w.f->getTo();

    // delete old edges
    _targetGraph->deleteEdge(w.e->getFrom(), w.e->getTo());
    _targetGraph->deleteEdge(w.f->getFrom(), w.f->getTo());

    assert(a != b);

    // add new edges
    _targetGraph->addNode(a);
    _targetGraph->addNode(b);
    Edge* eaE =_targetGraph->addEdge(wefrom, a);
    Edge* abE =_targetGraph->addEdge(a, b);
    Edge* ebE =_targetGraph->addEdge(b, weto);

    Edge* faE = 0;
    Edge* fbE = 0;

    if (reversed) {
      faE =_targetGraph->addEdge(a, wfto);
      fbE =_targetGraph->addEdge(wffrom, b);
    } else {
      faE =_targetGraph->addEdge(wffrom, a);
      fbE =_targetGraph->addEdge(b, wfto);
    }

    if (eaE) eaE->addEdgeTripGeom(eaEdgeGeom);
    if (abE) abE->addEdgeTripGeom(abEdgeGeom);
    if (ebE) ebE->addEdgeTripGeom(ecEdgeGeom);

    if (faE) faE->addEdgeTripGeom(faEdgeGeom);
    if (fbE) fbE->addEdgeTripGeom(fcEdgeGeom);
  }
}

// _____________________________________________________________________________
geo::PolyLine GraphBuilder::getSubPolyLine(Stop* a, Stop* b, Trip* t) {
  if (!t->getShape()) {
    return geo::PolyLine(getProjectedPoint(a->getLat(), a->getLng()),
      getProjectedPoint(b->getLat(), b->getLng()));
  }

  auto pl = _polyLines.find(t->getShape());
  if (pl == _polyLines.end()) {
    // generate polyline for this shape
    pl = _polyLines.insert(std::pair<gtfs::Shape*, geo::PolyLine>(t->getShape(), geo::PolyLine())).first;

    for (const auto& sp : t->getShape()->getPoints()) {
      pl->second << getProjectedPoint(sp.lat, sp.lng);
    }

    pl->second.simplify(10);
  }

  return pl->second.getSegment(getProjectedPoint(a->getLat(), a->getLng()),
    getProjectedPoint(b->getLat(), b->getLng()));
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

    n->setPos(util::geo::Point(x/c, y/c));
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

  util::geo::Point p = getProjectedPoint(curStop->getLat(), curStop->getLng());

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
         if (rg.getTripsUnordered().size() >= lc) {
           lc = rg.getTripsUnordered().size();
           g = &rg;
         }
      }

      if (g->getGeom().getLength() == 0) continue;

      NodeFront f(e, n, g);
      geo::PolyLine pl;

      if (g->getGeomDir() == n) {
        pl = g->getGeom().getOrthoLineAtDist(g->getGeom().getLength(), g->getTotalWidth());
      } else {
        pl = g->getGeom().getOrthoLineAtDist(0, g->getTotalWidth());
      }

      geo::PointOnLine curGeomPos = *(pl.getIntersections(f.refEtg->getGeom()).begin());
      if (g->getGeomDir() == n) {
        assert(curGeomPos.totalPos > .5);
      } else {
        assert(curGeomPos.totalPos < .5);
      }

      f.setGeom(pl);
      n->addMainDir(f);
    }

    // now, look at the nodes entire front geometries and expand them
    // until nothing overlaps
    //while (nodeHasOverlappingFronts(n)) {
    while (nodeHasOverlappingFronts(n)) {
      for (auto& f : n->getMainDirs()) {
        // TODO: store the reference ETG in the front
        geo::PointOnLine curGeomPos = *(f.geom.getIntersections(f.refEtg->getGeom()).begin());
        if (f.refEtg->getGeomDir() == n) {
          if (curGeomPos.totalPos < .5) goto exitloop;
          f.geom = f.refEtg->getGeom().getOrthoLineAtDist(curGeomPos.totalPos * f.refEtg->getGeom().getLength() - 5, f.refEtg->getTotalWidth());
        } else {
          if (curGeomPos.totalPos > .5) goto exitloop;
          f.geom = f.refEtg->getGeom().getOrthoLineAtDist(curGeomPos.totalPos * f.refEtg->getGeom().getLength() + 5, f.refEtg->getTotalWidth());
        }
      }
    }
    exitloop:
      continue;
  }
}

// _____________________________________________________________________________
bool GraphBuilder::nodeHasOverlappingFronts(const Node* n) const {
  for (size_t i = 0; i < n->getMainDirs().size(); ++i) {
    const NodeFront& fa = n->getMainDirs()[i];
    for (size_t j = 0; j < n->getMainDirs().size(); ++j) {
      const NodeFront& fb = n->getMainDirs()[j];
      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      if (n->getStops().size() > 0 && fa.geom.distTo(fb.geom) < (fa.refEtg->getSpacing() + fb.refEtg->getSpacing()) / 8) {
        return true;
      } else if (n->getStops().size() == 0 && fa.geom.distTo(fb.geom) < (fa.refEtg->getSpacing() + fb.refEtg->getSpacing()) * 2) {
        return true;
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
void GraphBuilder::freeNodes() {
  // TODO: move the creation of the node front geometries to the
  // method above. add method to "expand" node front geometries until they DO
  // NOT OVERLAP. after this, add an extension method to the polyline to double the
  // "stretch" the geometry to the desired length. use this as a cutting geom.
  for (auto n : *_targetGraph->getNodes()) {
    for (auto& f : n->getMainDirs()) {
      for (auto e : f.edges) {
        geo::PolyLine cutLine = f.geom;

        for (EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
          std::set<geo::PointOnLine, geo::PointOnLineCompare> iSects = cutLine.getIntersections(g.getGeom());
          if (iSects.size() > 0) {
            if (g.getGeomDir() !=n) {
              // cut at beginning
              g.setGeom(g.getGeom().getSegment(iSects.begin()->totalPos, 1));
            } else {
              // cut at end
              g.setGeom(g.getGeom().getSegment(0, (--iSects.end())->totalPos));
            }
          }
        }
      }
    }
  }
}
