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
using util::geo::SharedSegment;

namespace gtfs2topo {

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0){};
  ShrdSegWrap(Edge* e, Edge* f, SharedSegment s) : e(e), f(f), s(s){};
  Edge* e;
  Edge* f;
  SharedSegment s;
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

  bool lineDominatesSharedSeg(const ShrdSegWrap& w, Edge* e) const;

  // map of compiled polylines, to avoid calculating them each time
  std::unordered_map<gtfs::Shape*, PolyLine> _polyLines;

  Point getProjectedPoint(double lat, double lng, projPJ p) const;

  std::pair<bool, PolyLine> getSubPolyLine(const gtfs::Stop* a,
                                           const gtfs::Stop* b, gtfs::Trip* t,
                                           double distA, double distB);

  ShrdSegWrap getNextSharedSegment(BuildGraph* g, bool final, Grid<Edge*, Line>* grid) const;
  PolyLine getAveragedFromSharedSeg(const ShrdSegWrap& w) const;

  Node* addStop(const gtfs::Stop* curStop, uint8_t aggrLevel, BuildGraph* g);

  bool checkTripSanity(gtfs::Trip* t) const;
  bool checkShapeSanity(gtfs::Shape* t) const;

  bool combineNodes(Node* a, Node* b, BuildGraph* g);
  bool combineEdges(Edge* a, Edge* b, Node* n, BuildGraph* g);

  bool lineCrossesAtNode(const Node* a, const Edge* e, const Edge* f) const;

  Node* getNodeByStop(const BuildGraph* g, const gtfs::Stop* s, bool getParent) const;
  Node* getNodeByStop(const BuildGraph* g, const gtfs::Stop* s) const;
  Node* getNearestNode(const BuildGraph* g, const Point& p, double maxD) const;

  Box getGraphBoundingBox(const BuildGraph* g) const;
  Grid<Edge*, Line>  getGeoIndex(const BuildGraph* g) const;

  mutable std::set<const Edge*> _indEdges;
  mutable std::set<std::pair<const Edge*, const Edge*> > _indEdgesPairs;
  mutable std::map<std::pair<const Edge*, const Edge*>, size_t> _pEdges;
};

}  // namespace gtfs2topo

#endif  // GTFS2TOPO_BUILDER_BUILDER_H_
