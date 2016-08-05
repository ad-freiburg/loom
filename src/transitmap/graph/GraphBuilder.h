// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_GRAPHBUILDER_H_
#define TRANSITMAP_GRAPH_GRAPHBUILDER_H_

#include <algorithm>
#include <unordered_map>
#include <proj_api.h>
#include "transitgraph.h"
#include "gtfsparser/gtfs/feed.h"
#include "../geo/PolyLine.h"

namespace transitmapper {
namespace graph {

const static char* WGS84_PROJ = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0) {};
  ShrdSegWrap(Edge* e, Edge* f, geo::SharedSegment s) : e(e), f(f), s(s) {};
  Edge* e;
  Edge* f;
  geo::SharedSegment s;
};

class GraphBuilder {

 public:
  GraphBuilder(TransitGraph* targetGraph);

  void consume(const gtfs::Feed& f);
  void simplify();
  void createTopologicalNodes();
 private:
  TransitGraph* _targetGraph;
  projPJ _mercProj;

  // map of compiled polylines, to avoid calculating them each time
  std::unordered_map<gtfs::Shape*, geo::PolyLine> _polyLines;

  util::geo::Point getProjectedPoint(double lat, double lng) const;

  geo::PolyLine getSubPolyLine(gtfs::Stop* a, gtfs::Stop* b, gtfs::Trip* t);
  ShrdSegWrap getNextSharedSegment() const;
};

}}

#endif  // TRANSITMAP_GRAPH_GRAPHBUILDER_H_
