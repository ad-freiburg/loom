// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_BUILDER_BUILDER_H_
#define GTFS2TOPO_BUILDER_BUILDER_H_

#include <proj_api.h>
#include <algorithm>
#include <unordered_map>
#include "./../config/GraphBuilderConfig.h"
#include "./../graph/Graph.h"
#include "gtfsparser/gtfs/Feed.h"
#include "pbutil/geo/PolyLine.h"

using namespace gtfs2topo::graph;
using pbutil::geo::PolyLine;
using pbutil::geo::SharedSegment;

namespace gtfs2topo {

const static char* WGS84_PROJ =
    "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

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

  // build a graph from a gtfs feed
  void consume(const gtfs::Feed& f, Graph* g);

  // simpliyfy the graph
  void simplify(Graph* g);
  bool createTopologicalNodes(Graph* g, bool final);
  void averageNodePositions(Graph* g);
  void removeEdgeArtifacts(Graph* g);
  void removeNodeArtifacts(Graph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;

  bool lineDominatesSharedSeg(const ShrdSegWrap& w, Edge* e) const;

  // map of compiled polylines, to avoid calculating them each time
  std::unordered_map<gtfs::Shape*, PolyLine> _polyLines;

  Point getProjectedPoint(double lat, double lng, projPJ p) const;

  std::pair<bool, PolyLine> getSubPolyLine(gtfs::Stop* a, gtfs::Stop* b,
                                           gtfs::Trip* t, double distA,
                                           double distB, projPJ p);

  ShrdSegWrap getNextSharedSegment(Graph* g, bool final) const;
  PolyLine getAveragedFromSharedSeg(const ShrdSegWrap& w) const;

  Node* addStop(gtfs::Stop* curStop, uint8_t aggrLevel, Graph* g);

  bool checkTripSanity(gtfs::Trip* t) const;
  bool checkShapeSanity(gtfs::Shape* t) const;

  bool combineNodes(Node* a, Node* b, Graph* g);
  bool combineEdges(Edge* a, Edge* b, Node* n, Graph* g);

  bool lineCrossesAtNode(const Node* a, const Edge* e, const Edge* f) const;

  mutable std::set<const Edge*> _indEdges;
  mutable std::map<const Edge*, size_t> _pEdges;
};

}  // namespace gtfs2topo

#endif  // GTFS2TOPO_BUILDER_BUILDER_H_
