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
  // add all the nodes first. the TransitGraph maintains
  // a map stationid->nodeid for us

  // TODO: nicer access to internal iterators in feed, but
  // dont expose the map!
  // (maybe iterators to only the values???)
  for (auto s = f.stopsBegin(); s != f.stopsEnd(); ++s) {
    Stop* curStop = s->second;

    // reproject to graph projection
    double x = curStop->getLng();
    double y = curStop->getLat();
    x *= DEG_TO_RAD;
    y *= DEG_TO_RAD;

    pj_transform(_mercProj, _targetGraph->getProjection(), 1, 1, &x, &y, 0);

    _targetGraph->addNode(new Node(x, y, curStop));
  }
}
