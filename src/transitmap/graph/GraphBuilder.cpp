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
    //if (t->second->getRoute()->getShortName() != "3") continue;
    //if (!(t->second->getShape())) continue;
    auto st = t->second->getStopTimes().begin();

    StopTime prev = *st;
    addStop(prev.getStop(), AGGREGATE_STOPS);
    ++st;

    for (; st != t->second->getStopTimes().end(); ++st) {
      const StopTime& cur = *st;
      addStop(cur.getStop(), AGGREGATE_STOPS);

      Node* fromNode = _targetGraph->getNodeByStop(
        prev.getStop(),
        AGGREGATE_STOPS
      );
      Node* toNode =  _targetGraph->getNodeByStop(
        cur.getStop(),
        AGGREGATE_STOPS
      );

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
      if (_indEdges.find(e) != _indEdges.end() || e->getEdgeTripGeoms()->size() != 1) continue;
      // TODO: outfactor this _______
      for (auto nt : *_targetGraph->getNodes()) {
        for (auto toTest : nt->getAdjListOut()) {
          // TODO: only check edges with a SINGLE geometry atm, see also check above
          if (_indEdges.find(toTest) != _indEdges.end() || toTest->getEdgeTripGeoms()->size() != 1) continue;
          if (e != toTest) {
            geo::SharedSegments s = e->getEdgeTripGeoms()->front().getGeom().getSharedSegments(toTest->getEdgeTripGeoms()->front().getGeom());
            i++;
            if (s.segments.size() > 0 && util::geo::dist(s.segments[0].first.p, s.segments[0].second.p) > 50) {
              _pEdges[e]++;
              if (_pEdges[e] > 100) {
                std::cout << "Too many optimiziations for " << e << ", preventing further..." << std::endl;
                _indEdges.insert(e);
              }
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
    if (fap.totalPos > fbp.totalPos) {
      geo::PointOnLine temp = fap;
      fap = fbp;
      fbp = temp;
      reversed = true;
    }

    geo::PolyLine fa = compEdgeGeom.getGeom().getSegment(0, fap.totalPos);
    geo::PolyLine fab = compEdgeGeom.getGeom().getSegment(fap, fbp);
    geo::PolyLine fc = compEdgeGeom.getGeom().getSegment(fbp.totalPos, 1);

    if (reversed) {
      fab.reverse();
    }
    std::vector<const geo::PolyLine*> avg;
    avg.push_back(&ab);
    avg.push_back(&fab);
    ab = geo::PolyLine::average(avg);

    // new node at the start of the shared segment
    Node* a = 0;

    if (util::geo::dist(curEdgeGeom.getGeom().projectOn(w.e->getFrom()->getPos()).p, w.s.first.p) < 50) {
      a = w.e->getFrom();
    } else if (util::geo::dist(curEdgeGeom.getGeom().projectOn(w.e->getTo()->getPos()).p, w.s.first.p) < 50) {
      a = w.e->getTo();
    } else {
      a = new Node(w.s.first.p);
    }

    // new node at the end of the shared segment
    Node* b = 0;

    if (util::geo::dist(curEdgeGeom.getGeom().projectOn(w.e->getTo()->getPos()).p, w.s.second.p) < 50) {
      b = w.e->getTo();
    } else if (util::geo::dist(curEdgeGeom.getGeom().projectOn(w.e->getFrom()->getPos()).p, w.s.second.p) < 50) {
      b = w.e->getFrom();
    } else {
      b = new Node(w.s.second.p);
    }

    EdgeTripGeom eaEdgeGeom(ea, a);
    EdgeTripGeom abEdgeGeom(ab, b);
    EdgeTripGeom ecEdgeGeom(ec, w.e->getTo());

    EdgeTripGeom faEdgeGeom(fa, a);
    EdgeTripGeom fcEdgeGeom(fc, w.f->getTo());

    for (auto& r : curEdgeGeom.getTrips()) {
      for (auto& t : r.second.trips) {
        if (r.second.direction == w.e->getTo()) {
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

    for (auto& r : compEdgeGeom.getTrips()) {
      for (auto& t : r.second.trips) {
        if (r.second.direction == w.f->getTo()) {
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


    // add new edges
    _targetGraph->addNode(a);
    _targetGraph->addNode(b);
    Edge* eaE =_targetGraph->addEdge(w.e->getFrom(), a);
    Edge* abE =_targetGraph->addEdge(a, b);
    Edge* ebE =_targetGraph->addEdge(b, w.e->getTo());

    Edge* faE = 0;
    Edge* fbE = 0;

    if (reversed) {
      faE =_targetGraph->addEdge(w.f->getFrom(), b);
      fbE =_targetGraph->addEdge(a, w.f->getTo());
    } else {
      faE =_targetGraph->addEdge(w.f->getFrom(), a);
      fbE =_targetGraph->addEdge(b, w.f->getTo());
    }

    if (eaE) eaE->addEdgeTripGeom(eaEdgeGeom);
    if (abE) abE->addEdgeTripGeom(abEdgeGeom);
    if (ebE) ebE->addEdgeTripGeom(ecEdgeGeom);

    if (faE) faE->addEdgeTripGeom(faEdgeGeom);
    if (fbE) fbE->addEdgeTripGeom(fcEdgeGeom);

    // delete old edges
    _targetGraph->deleteEdge(w.e->getFrom(), w.e->getTo());
    _targetGraph->deleteEdge(w.f->getFrom(), w.f->getTo());
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

    pl->second.simplify(3);
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
        if (eg.getGeomDir() == e->getTo()) {
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
        if (eg.getGeomDir() == e->getTo()) {
          x += eg.getGeom().getLine().back().get<0>();
          y += eg.getGeom().getLine().back().get<1>();
        } else {
          x += eg.getGeom().getLine().front().get<0>();
          y += eg.getGeom().getLine().front().get<1>();
        }
        c++;
      }
    }

    n->setPos(util::geo::Point(x/c, y/c));
  }
}

// _____________________________________________________________________________
void GraphBuilder::addStop(gtfs::Stop* curStop, uint8_t aggrLevel) {
  if (aggrLevel && curStop->getParentStation() != 0) {
    addStop(curStop->getParentStation(), aggrLevel);
  }

  util::geo::Point p = getProjectedPoint(curStop->getLat(), curStop->getLng());

  Node* n = 0;

  if (aggrLevel > 1) {
    n = _targetGraph->getNearestNode(p, 100);
  }

  if (n) {
    n->addStop(curStop);
  } else {
    _targetGraph->addNode(
      new Node(
        p,
        curStop
      )
    );
  }
}

// _____________________________________________________________________________
void GraphBuilder::writeMainDirs() {
  for (auto n : *_targetGraph->getNodes()) {
    std::set<Edge*> eSet;
    eSet.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
    eSet.insert(n->getAdjListOut().begin(), n->getAdjListOut().end());

    bool first = true;
    int32_t ref;

    for (Edge* e : eSet) {
      if (e->getEdgeTripGeoms()->size() == 0) continue;

      // atm, always take the first edge trip geometry
      const EdgeTripGeom& g = e->getEdgeTripGeoms()->front();
      int32_t angle;
      util::geo::Point p = g.getGeom().projectOn(n->getPos()).p;
      if (g.getGeomDir() == n) {
        angle = util::geo::angBetween(p, g.getGeom().getPointAtDist(g.getGeom().getLength() - 50).p);
      } else {
        angle = util::geo::angBetween(p, g.getGeom().getPointAtDist(50).p);
      }

      /**

      angle = ((angle % 360) + 360) % 360;

      if (first) {
        first = false;
        ref = angle;
      }

      int step = 90;
      div_t res;

      int32_t diff = (((angle - ref) % 360) + 360) % 360;
      res = div(diff, step);

      if (res.rem < step/3) angle -= res.rem;
      if (res.rem >= 2*(step/3)) angle += step - res.rem;
      **/

      n->addMainDir(NodeFront((angle + 90) / (180 / M_PI), e));
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::freeNodes(double d) {
  for (auto n : *_targetGraph->getNodes()) {
    size_t c = 0;
    for (auto f : n->getMainDirs()) {
      size_t curN = 0;
      for (auto e : f.edges) {
        for (auto g : *e->getEdgeTripGeoms()) {
          curN += g.getTrips()->size();
        }
      }
      if (curN > c) c = curN;
    }
    std::cout << c << std::endl;
    for (auto f : n->getMainDirs()) {
      for (auto e : f.edges) {
        for (auto g : *e->getEdgeTripGeoms()) {
          if (g.getGeomDir() != n) {
            g.setGeom(g.getGeom().getSegment(g.getGeom().getPointAtDist(d * c).totalPos, 1));
          } else {
            g.setGeom(g.getGeom().getSegment(0, g.getGeom().getPointAtDist(g.getGeom().getLength() - d * c).totalPos));
          }
        }
      }
    }
  }
}
