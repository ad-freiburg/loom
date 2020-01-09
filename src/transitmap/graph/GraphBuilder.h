// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_GRAPHBUILDER_H_
#define TRANSITMAP_GRAPH_GRAPHBUILDER_H_

#include <proj_api.h>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <vector>
#include "TransitGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "util/geo/PolyLine.h"

namespace transitmapper {
namespace graph {

using util::geo::SharedSegment;

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0){};
  ShrdSegWrap(shared::linegraph::LineEdge* e, shared::linegraph::LineEdge* f,
              SharedSegment<double> s)
      : e(e), f(f), s(s){};
  shared::linegraph::LineEdge* e;
  shared::linegraph::LineEdge* f;
  SharedSegment<double> s;
};

class GraphBuilder {
 public:
  GraphBuilder(const config::Config* cfg);
  bool build(std::istream* s, graph::TransitGraph* g);

  void writeMainDirs(TransitGraph* g);
  void expandOverlappinFronts(TransitGraph* g);
  void writeInitialConfig(TransitGraph* g);
  void createMetaNodes(TransitGraph* g);

 private:
  const config::Config* _cfg;
  projPJ _mercProj;

  std::set<shared::linegraph::NodeFront*> nodeGetOverlappingFronts(const TransitGraph* g,
      const shared::linegraph::LineNode* n) const;
  void freeNodeFront(const shared::linegraph::LineNode* n, shared::linegraph::NodeFront* f);

  std::vector<shared::linegraph::NodeFront> getNextMetaNodeCand(
      TransitGraph* g) const;
  std::vector<shared::linegraph::NodeFront> getOpenNodeFronts(
      const graph::TransitGraph* g, const shared::linegraph::LineNode* n) const;
  std::vector<shared::linegraph::NodeFront> getClosedNodeFronts(
      const graph::TransitGraph* g, const shared::linegraph::LineNode* n) const;
  bool isClique(const graph::TransitGraph* g, std::set<const shared::linegraph::LineNode*> potClique) const;

  bool nodeFrontsOverlap(const TransitGraph* g, const shared::linegraph::NodeFront& a,
                         const shared::linegraph::NodeFront& b) const;
  mutable std::set<const shared::linegraph::LineEdge*> _indEdges;
  mutable std::map<const shared::linegraph::LineEdge*, size_t> _pEdges;
};

}  // namespace graph
}  // namespace transitmapper

#endif  // TRANSITMAP_GRAPH_GRAPHBUILDER_H_
