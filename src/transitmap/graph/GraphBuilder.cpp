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

  bool AGGREGATE_STOPS = true;

  // add all the nodes first. the TransitGraph maintains
  // a map stationid->nodeid for us

  // TODO: nicer access to internal iterators in feed, but
  // dont expose the map!
  for (auto s = f.stopsBegin(); s != f.stopsEnd(); ++s) {
    Stop* curStop = s->second;
    if (AGGREGATE_STOPS && curStop->getParentStation() != 0) continue;

    // reproject to graph projection
    double x = curStop->getLng();
    double y = curStop->getLat();
    x *= DEG_TO_RAD;
    y *= DEG_TO_RAD;

    pj_transform(_mercProj, _targetGraph->getProjection(), 1, 1, &x, &y, 0);

    _targetGraph->addNode(
      new Node(
        getProjectedPoint(curStop->getLat(), curStop->getLng()),
        curStop
      )
    );
  }

  size_t cur = 0;
  for (auto t = f.tripsBegin(); t != f.tripsEnd(); ++t) {
    cur++;
    if (t->second->getStopTimes().size() < 2) continue;

    if (t->second->getRoute()->getType() != gtfs::Route::TYPE::TRAM) continue;
    if (!(t->second->getShape())) continue;

    std::cout << "AT " << cur << std::endl;

    auto st = t->second->getStopTimes().begin();

    StopTime prev = *st;
    ++st;
    for (; st != t->second->getStopTimes().end(); ++st) {
      const StopTime& cur = *st;
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
