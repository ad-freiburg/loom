// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_GRAPH_H_
#define OCTI_COMBGRAPH_GRAPH_H_

#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"
#include "shared/linegraph/LineGraph.h"
#include "util/graph/UndirGraph.h"

namespace octi {
namespace combgraph {

typedef util::graph::Node<CombNodePL, CombEdgePL> CombNode;
typedef util::graph::Edge<CombNodePL, CombEdgePL> CombEdge;

using shared::linegraph::LineGraph;
using octi::combgraph::EdgeOrdering;

class CombGraph : public util::graph::UndirGraph<CombNodePL, CombEdgePL> {
 public:
  CombGraph(const LineGraph* g);
  CombGraph(const LineGraph* g, bool collapse);

  EdgeOrdering getEdgeOrderingForNode(CombNode* n) const;
  EdgeOrdering getEdgeOrderingForNode(
      CombNode* n, bool useOrigNextNode,
      const std::map<CombNode*, util::geo::DPoint>& newPos) const;

 private:
  void build(const LineGraph* source);
  void combineDeg2();
  void writeEdgeOrdering();
};

}  // combgraph
}  // octi

#endif  // OCTI_COMBGRAPH_GRAPH_H_
