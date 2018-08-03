// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_BUILDER_BUILDER_H_
#define GTFS2TOPO_BUILDER_BUILDER_H_

#include <proj_api.h>
#include <algorithm>
#include <unordered_map>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "gtfs2topo/config/GraphBuilderConfig.h"
#include "util/graph/Graph.h"
#include "util/geo/PolyLine.h"
#include "util/geo/Grid.h"
#include "util/geo/Geo.h"
#include "gtfs2topo/graph/BuildGraph.h"

using namespace gtfs2topo::graph;
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

namespace gtfs2topo {

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0){};
  ShrdSegWrap(Edge* e, Edge* f, SharedSegment<double> s) : e(e), f(f), s(s){};
  Edge* e;
  Edge* f;
  SharedSegment<double> s;
};

class Builder {
 public:
  Builder(const config::Config* cfg);

  // build a BuildGraph from a gtfs feed
  void consume(const gtfs::Feed& f, BuildGraph* g);

  // simpliyfy the BuildGraph
  void simplify(BuildGraph* g);
  bool createTopologicalNodes(BuildGraph* g, bool final);
  void averageNodePositions(BuildGraph* g);
  void removeEdgeArtifacts(BuildGraph* g);
  void removeNodeArtifacts(BuildGraph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;
  projPJ _graphProj;

  std::map<const gtfs::Stop*, Node*> _stopNodes;

  bool lineDominatesSharedSeg(const ShrdSegWrap& w, Edge* e) const;

  // map of compiled polylines, to avoid calculating them each time
  std::unordered_map<gtfs::Shape*, PolyLine<double>> _polyLines;

  DPoint getProjectedPoint(double lat, double lng, projPJ p) const;

  std::pair<bool, PolyLine<double>> getSubPolyLine(const gtfs::Stop* a,
                                           const gtfs::Stop* b, gtfs::Trip* t,
                                           double distA, double distB);

  ShrdSegWrap getNextSharedSegment(BuildGraph* g, bool final, EdgeGrid* grid) const;
  PolyLine<double> getAveragedFromSharedSeg(const ShrdSegWrap& w) const;

  Node* addStop(const gtfs::Stop* curStop, uint8_t aggrLevel, BuildGraph* g, NodeGrid* grid);

  bool checkTripSanity(gtfs::Trip* t) const;
  bool checkShapeSanity(gtfs::Shape* t) const;

  bool combineNodes(Node* a, Node* b, BuildGraph* g);
  bool combineEdges(Edge* a, Edge* b, Node* n, BuildGraph* g);

  bool lineCrossesAtNode(const Node* a, const Edge* e, const Edge* f) const;

  Node* getNodeByStop(const BuildGraph* g, const gtfs::Stop* s, bool getParent) const;
  Node* getNodeByStop(const BuildGraph* g, const gtfs::Stop* s) const;
  Node* getNearestStop(const BuildGraph* g, const DPoint& p, double maxD, const NodeGrid* grid) const;

  DBox getGraphBoundingBox(const BuildGraph* g) const;
  EdgeGrid getGeoIndex(const BuildGraph* g) const;

  mutable std::set<const Edge*> _indEdges;
  mutable std::set<std::pair<const Edge*, const Edge*> > _indEdgesPairs;
  mutable std::map<std::pair<const Edge*, const Edge*>, size_t> _pEdges;
};

}  // namespace gtfs2topo

#endif  // GTFS2TOPO_BUILDER_BUILDER_H_
