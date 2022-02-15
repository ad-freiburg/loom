// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_GRAPHBUILDER_H_
#define TRANSITMAP_GRAPH_GRAPHBUILDER_H_

#include <algorithm>
#include <set>
#include <unordered_map>
#include <vector>
#include "shared/rendergraph/RenderGraph.h"
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

  void writeNodeFronts(shared::rendergraph::RenderGraph* g);
  void expandOverlappinFronts(shared::rendergraph::RenderGraph* g);

 private:
  const config::Config* _cfg;

  std::set<shared::linegraph::NodeFront*> nodeGetOverlappingFronts(
      const shared::rendergraph::RenderGraph* g,
      const shared::linegraph::LineNode* n) const;
  void freeNodeFront(const shared::linegraph::LineNode* n,
                     shared::linegraph::NodeFront* f);

  bool nodeFrontsOverlap(const shared::rendergraph::RenderGraph* g,
                         const shared::linegraph::NodeFront& a,
                         const shared::linegraph::NodeFront& b,
                         double d) const;
  mutable std::set<const shared::linegraph::LineEdge*> _indEdges;
  mutable std::map<const shared::linegraph::LineEdge*, size_t> _pEdges;
};

}  // namespace graph
}  // namespace transitmapper

#endif  // TRANSITMAP_GRAPH_GRAPHBUILDER_H_
