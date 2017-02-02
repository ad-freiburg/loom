// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SKELETONBUILDER_BUILDER_BUILDER_H_
#define SKELETONBUILDER_BUILDER_BUILDER_H_

#include <algorithm>
#include <unordered_map>
#include <proj_api.h>
#include "transitmap/graph/TransitGraph.h"
#include "gtfsparser/gtfs/Feed.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/geo/PolyLine.h"

using namespace transitmapper;
using namespace graph;

namespace skeletonbuilder {

const static char* WGS84_PROJ = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0) {};
  ShrdSegWrap(Edge* e, Edge* f, geo::SharedSegment s) : e(e), f(f), s(s) {};
  Edge* e;
  Edge* f;
  geo::SharedSegment s;
};

class Builder {
 public:
  Builder(const config::Config* cfg);

  // build a graph from a gtfs feed
  void consume(const gtfs::Feed& f, TransitGraph* targetGraph);

  // simpliyfy the graph
  void simplify(TransitGraph* g);
  bool createTopologicalNodes(TransitGraph* g);
  void averageNodePositions(TransitGraph* g);
  void removeArtifacts(TransitGraph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;

  bool lineDominatesSharedSeg(const ShrdSegWrap& w, Edge* e) const;

  // map of compiled polylines, to avoid calculating them each time
  std::unordered_map<gtfs::Shape*, geo::PolyLine> _polyLines;

  util::geo::Point getProjectedPoint(double lat, double lng, projPJ p) const;

  std::pair<bool, geo::PolyLine> getSubPolyLine(gtfs::Stop* a, gtfs::Stop* b,
      gtfs::Trip* t, double distA, double distB, projPJ p);

  ShrdSegWrap getNextSharedSegment(TransitGraph* g) const;
  geo::PolyLine getAveragedFromSharedSeg(const ShrdSegWrap& w) const;

  std::set<NodeFront*> nodeGetOverlappingFronts(const Node* n) const;

  Node* addStop(gtfs::Stop* curStop, uint8_t aggrLevel, TransitGraph* g);
  void freeNodeFront(NodeFront* f);

  bool nodeFrontsOverlap(const NodeFront& a, const NodeFront& b) const;

  bool checkTripSanity(gtfs::Trip* t) const;
  bool checkShapeSanity(gtfs::Shape* t) const;

  void combineNodes(Node* a, Node* b, TransitGraph* g);

  mutable std::set<const Edge*> _indEdges;
  mutable std::map<const Edge*, size_t> _pEdges;
};

}  // namespace skeletonbuilder

#endif  // SKELETONBUILDER_BUILDER_BUILDER_H_

