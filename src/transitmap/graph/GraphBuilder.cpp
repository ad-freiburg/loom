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
  // (maybe iterators to only the values???)
  for (auto s = f.stopsBegin(); s != f.stopsEnd(); ++s) {
    Stop* curStop = s->second;
    if (AGGREGATE_STOPS && curStop->getParentStation() != 0) continue;

    // reproject to graph projection
    double x = curStop->getLng();
    double y = curStop->getLat();
    x *= DEG_TO_RAD;
    y *= DEG_TO_RAD;

    pj_transform(_mercProj, _targetGraph->getProjection(), 1, 1, &x, &y, 0);

    _targetGraph->addNode(new Node(x, y, curStop));
  }

  for (auto t = f.tripsBegin(); t != f.tripsEnd(); ++t) {
    if (t->second->getStopTimes().size() < 2) continue;

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

      exE->addTrip(t->second, geo::PolyLine(fromNode->getPos(), toNode->getPos()));
      prev = cur;
    }
  }
}

// _____________________________________________________________________________ 
void GraphBuilder::simplify() {
  // try to merge both-direction edges into a single one

  for (auto n : *_targetGraph->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      // check if counter-edge exists
      Edge* ce = targetGraph->getEdge(e->getTo(), e->getFrom());
      
      for (auto e 
    }
  }
}
