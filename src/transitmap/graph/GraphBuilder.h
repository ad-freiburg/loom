// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_GRAPHBUILDER_H_
#define TRANSITMAP_GRAPH_GRAPHBUILDER_H_

#include <algorithm>
#include <unordered_map>
#include <proj_api.h>
#include "TransitGraph.h"
#include "gtfsparser/gtfs/Feed.h"
#include "./../config/TransitMapConfig.h"
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
  GraphBuilder(const config::Config* cfg);
  bool build(std::istream* s, graph::TransitGraph* g);

  void writeMainDirs(TransitGraph* g);
  void expandOverlappinFronts(TransitGraph* g);
  void writeInitialConfig(TransitGraph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;

  std::set<NodeFront*> nodeGetOverlappingFronts(const Node* n) const;
  void freeNodeFront(NodeFront* f);

  bool nodeFrontsOverlap(const NodeFront& a, const NodeFront& b) const;
  mutable std::set<const Edge*> _indEdges;
  mutable std::map<const Edge*, size_t> _pEdges;
};

}  // namespace graph
}  // namespace transitmapper

#endif  // TRANSITMAP_GRAPH_GRAPHBUILDER_H_
