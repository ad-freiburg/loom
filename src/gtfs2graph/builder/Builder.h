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

using namespace gtfs2graph::graph;
using namespace ad::cppgtfs;
using util::geo::Grid;
using util::geo::Box;
using util::geo::Line;
using util::geo::PolyLine;
using util::geo::Point;
using util::geo::Line;
using util::geo::DPoint;
using util::geo::DLine;
using util::geo::DBox;
using util::geo::SharedSegment;

typedef Grid<Node*, Point, double> NodeGrid;
typedef Grid<Edge*, Line, double> EdgeGrid;

namespace gtfs2graph {

class Builder {
 public:
  Builder(const config::Config* cfg);

  // build a BuildGraph from a gtfs feed
  void consume(const gtfs::Feed& f, BuildGraph* g);

  // simpliyfy the BuildGraph
  void simplify(BuildGraph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;
  projPJ _graphProj;

  std::map<const gtfs::Stop*, Node*> _stopNodes;

  // map of compiled polylines, to avoid calculating them each time
  std::unordered_map<gtfs::Shape*, PolyLine<double>> _polyLines;

  DPoint getProjectedPoint(double lat, double lng, projPJ p) const;

  std::pair<bool, PolyLine<double>> getSubPolyLine(const gtfs::Stop* a,
                                                   const gtfs::Stop* b,
                                                   gtfs::Trip* t, double distA,
                                                   double distB);

  Node* addStop(const gtfs::Stop* curStop, uint8_t aggrLevel, BuildGraph* g,
                NodeGrid* grid);

  bool checkTripSanity(gtfs::Trip* t) const;
  bool checkShapeSanity(gtfs::Shape* t) const;

  Node* getNodeByStop(const BuildGraph* g, const gtfs::Stop* s,
                      bool getParent) const;
  Node* getNodeByStop(const BuildGraph* g, const gtfs::Stop* s) const;
  Node* getNearestStop(const DPoint& p, double maxD,
                       const NodeGrid* grid) const;

  mutable std::set<const Edge*> _indEdges;
  mutable std::set<std::pair<const Edge*, const Edge*>> _indEdgesPairs;
  mutable std::map<std::pair<const Edge*, const Edge*>, size_t> _pEdges;
};

}  // namespace gtfs2graph

#endif  // GTFS2GRAPH_BUILDER_BUILDER_H_
