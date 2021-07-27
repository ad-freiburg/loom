// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_GRAPH_GRAPHBUILDER_H_
#define LOOMP_GRAPH_GRAPHBUILDER_H_

#include <algorithm>
#include <set>
#include <unordered_map>
#include <vector>
#include "loom/config/TransitMapConfig.h"
#include "shared/rendergraph/RenderGraph.h"
#include "util/geo/PolyLine.h"

namespace loom {
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

  void writeMainDirs(shared::rendergraph::RenderGraph* g);
  void writeInitialConfig(shared::rendergraph::RenderGraph* g);

 private:
  const config::Config* _cfg;

  std::set<shared::linegraph::NodeFront*> nodeGetOverlappingFronts(
      const shared::rendergraph::RenderGraph* g,
      const shared::linegraph::LineNode* n) const;
  void freeNodeFront(const shared::linegraph::LineNode* n,
                     shared::linegraph::NodeFront* f);

  bool nodeFrontsOverlap(const shared::rendergraph::RenderGraph* g,
                         const shared::linegraph::NodeFront& a,
                         const shared::linegraph::NodeFront& b) const;
  mutable std::set<const shared::linegraph::LineEdge*> _indEdges;
  mutable std::map<const shared::linegraph::LineEdge*, size_t> _pEdges;
};

}  // namespace graph
}  // namespace loom

#endif  // LOOM_GRAPH_GRAPHBUILDER_H_
