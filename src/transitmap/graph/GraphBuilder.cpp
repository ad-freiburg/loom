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

  // add all the nodes first. the TransitGraph maintains
  // a map stationid->nodeid for us

  // TODO: nicer access to internal iterators in feed, but
  // dont expose the map!
  for (auto s = f.stopsBegin(); s != f.stopsEnd(); ++s) {
    Stop* curStop = s->second;
    if (AGGREGATE_STOPS && curStop->getParentStation() != 0) continue;

    util::geo::Point p = getProjectedPoint(curStop->getLat(), curStop->getLng());

    std::cout << boost::geometry::wkt(p) << std::endl;
    Node* n = 0;

    if (AGGREGATE_STOPS > 1) {
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

  size_t cur = 0;
  for (auto t = f.tripsBegin(); t != f.tripsEnd(); ++t) {
    cur++;
    if (t->second->getStopTimes().size() < 2) continue;

    if (t->second->getRoute()->getType() == gtfs::Route::TYPE::BUS) continue;
    //if (t->second->getRoute()->getShortName() != "3") continue;
    //if (!(t->second->getShape())) continue;
    auto st = t->second->getStopTimes().begin();

    StopTime prev = *st;
    ++st;
    /**
      while ((!prev.getStop()->getParentStation() || (prev.getStop()->getParentStation()->getId() != "Parent30501" && 
          prev.getStop()->getParentStation()->getId() != "Parent30505")) && st != t->second->getStopTimes().end()) {
        prev = *st;
        ++st;
      }
**/
    for (; st != t->second->getStopTimes().end(); ++st) {
      const StopTime& cur = *st;
      /**
      if (!cur.getStop()->getParentStation() || (cur.getStop()->getParentStation()->getId() != "Parent30501" && 
          cur.getStop()->getParentStation()->getId() != "Parent30505")) continue;
      **/
      Node* fromNode = _targetGraph->getNodeByStop(
        prev.getStop(),
        AGGREGATE_STOPS
      );
      Node* toNode =  _targetGraph->getNodeByStop(
        cur.getStop(),
        AGGREGATE_STOPS
      );

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
  for (auto n : *_targetGraph->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      if (e->getEdgeTripGeoms()->size() != 1) continue;
      // TODO: outfactor this _______
      for (auto nt : *_targetGraph->getNodes()) {
        for (auto toTest : nt->getAdjListOut()) {
          // TODO: only check edges with a SINGLE geometry atm, see also check above
          if (toTest->getEdgeTripGeoms()->size() != 1) continue;
          if (e != toTest) {
            geo::SharedSegments s = e->getEdgeTripGeoms()->front().getGeom().getSharedSegments(toTest->getEdgeTripGeoms()->front().getGeom());
            if (s.segments.size() > 0 && boost::geometry::distance(s.segments.front().first.p, s.segments.front().second.p) > 10) {
              return ShrdSegWrap(e, toTest, s.segments.front());
            }
          }
        }
      }
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

    if (reversed) fab.reverse();

    std::vector<const geo::PolyLine*> avg;
    avg.push_back(&ab);
    avg.push_back(&fab);
    // ab = geo::PolyLine::average(avg);

    // new node at the start of the shared segment
    Node* a = new Node(w.s.first.p);

    // new node at the end of the shared segment
    Node* b = new Node(w.s.second.p);


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

    Edge* faE =_targetGraph->addEdge(w.f->getFrom(), a);
    Edge* fbE =_targetGraph->addEdge(b, w.f->getTo());

    eaE->addEdgeTripGeom(eaEdgeGeom);
    abE->addEdgeTripGeom(abEdgeGeom);
    ebE->addEdgeTripGeom(ecEdgeGeom);

    faE->addEdgeTripGeom(faEdgeGeom);
    fbE->addEdgeTripGeom(fcEdgeGeom);

    std::cout << w.e << std::endl;

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
  }

  return pl->second.getSegment(getProjectedPoint(a->getLat(), a->getLng()),
    getProjectedPoint(b->getLat(), b->getLng()));
}
