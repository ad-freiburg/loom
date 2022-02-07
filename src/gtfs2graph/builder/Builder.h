// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2GRAPH_BUILDER_BUILDER_H_
#define GTFS2GRAPH_BUILDER_BUILDER_H_

#include <proj_api.h>
#include <algorithm>
#include <unordered_map>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "gtfs2graph/config/GraphBuilderConfig.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Graph.h"

using gtfs2graph::graph::BuildGraph;
using gtfs2graph::graph::Edge;
using gtfs2graph::graph::Node;
using util::geo::Box;
using util::geo::DBox;
using util::geo::DLine;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Line;
using util::geo::Point;
using util::geo::PolyLine;
using util::geo::SharedSegment;

typedef Grid<Node*, Point, double> NodeGrid;
typedef Grid<Edge*, Line, double> EdgeGrid;

namespace gtfs2graph {

class Builder {
 public:
  Builder(const config::Config* cfg);

  // build a BuildGraph from a gtfs feed
  void consume(const ad::cppgtfs::gtfs::Feed& f, BuildGraph* g);

  // simplify the BuildGraph
  void simplify(BuildGraph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;
  projPJ _graphProj;

  std::map<const ad::cppgtfs::gtfs::Stop*, Node*> _stopNodes;

  // map of compiled polylines, to avoid calculating them each time
  std::unordered_map<ad::cppgtfs::gtfs::Shape*, PolyLine<double>> _polyLines;

  DPoint getProjectedPoint(double lat, double lng, projPJ p) const;

  std::pair<bool, PolyLine<double>> getSubPolyLine(
      const ad::cppgtfs::gtfs::Stop* a, const ad::cppgtfs::gtfs::Stop* b,
      ad::cppgtfs::gtfs::Trip* t, double distA, double distB);

  Node* addStop(const ad::cppgtfs::gtfs::Stop* curStop, BuildGraph* g,
                NodeGrid* grid);

  Node* getNodeByStop(const BuildGraph* g,
                      const ad::cppgtfs::gtfs::Stop* s) const;
};

}  // namespace gtfs2graph

#endif  // GTFS2GRAPH_BUILDER_BUILDER_H_
