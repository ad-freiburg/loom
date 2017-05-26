// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_GRAPHBUILDER_H_
#define TRANSITMAP_GRAPH_GRAPHBUILDER_H_

#include <algorithm>
#include <unordered_map>
#include <vector>
#include <set>
#include <proj_api.h>
#include "pbutil/geo/PolyLine.h"
#include "TransitGraph.h"
#include "./../config/TransitMapConfig.h"

namespace transitmapper {
namespace graph {

using pbutil::geo::SharedSegment;

const static char* WGS84_PROJ = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0) {};
  ShrdSegWrap(Edge* e, Edge* f, SharedSegment s) : e(e), f(f), s(s) {};
  Edge* e;
  Edge* f;
  SharedSegment s;
};

class GraphBuilder {
 public:
  GraphBuilder(const config::Config* cfg);
  bool build(std::istream* s, graph::TransitGraph* g);

  void combinePartnerRoutes(graph::TransitGraph* g);

  void writeMainDirs(TransitGraph* g);
  void expandOverlappinFronts(TransitGraph* g);
  void writeInitialConfig(TransitGraph* g);
  void createMetaNodes(TransitGraph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;

  std::set<NodeFront*> nodeGetOverlappingFronts(const Node* n) const;
  void freeNodeFront(NodeFront* f);

  std::map<const Route*, std::set<const Route*> >
  getPartnerRoutes(graph::TransitGraph* g) const;

  std::vector<NodeFront> getNextMetaNodeCand(TransitGraph* g) const;
  std::vector<NodeFront> getOpenNodeFronts(const Node* n) const;
  std::vector<NodeFront> getClosedNodeFronts(const Node* n) const;
  bool isClique(std::set<const Node*> potClique) const;

  bool nodeFrontsOverlap(const NodeFront& a, const NodeFront& b) const;
  mutable std::set<const Edge*> _indEdges;
  mutable std::map<const Edge*, size_t> _pEdges;
};

}  // namespace graph
}  // namespace transitmapper

#endif  // TRANSITMAP_GRAPH_GRAPHBUILDER_H_
