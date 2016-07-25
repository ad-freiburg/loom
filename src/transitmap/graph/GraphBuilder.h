// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_GRAPHBUILDER_H_
#define TRANSITMAP_GRAPH_GRAPHBUILDER_H_

#include <proj_api.h>
#include "transitgraph.h"
#include "gtfsparser/gtfs/feed.h"

namespace transitmapper {
namespace graph {

const static char* WGS84_PROJ = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

class GraphBuilder {

 public:
  GraphBuilder(TransitGraph* targetGraph);

  void consume(const gtfs::Feed& f);
  void simplify();
 private:
  TransitGraph* _targetGraph;
  projPJ _mercProj;
};

}}

#endif  // TRANSITMAP_GRAPH_GRAPHBUILDER_H_
