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

using octi::combgraph::EdgeOrdering;
using shared::linegraph::LineGraph;

class CombGraph : public util::graph::UndirGraph<CombNodePL, CombEdgePL> {
 public:
  CombGraph(const LineGraph* g);
  CombGraph(const LineGraph* g, bool collapse);

  EdgeOrdering getEdgeOrderingForNode(CombNode* n) const;
  EdgeOrdering getEdgeOrderingForNode(CombNode* n, bool useOrigNextNode) const;

  const util::geo::DBox& getBBox() const;

 private:
  util::geo::Box<double> _bbox;
  void build(const LineGraph* source);
  void combineDeg2();
  void writeEdgeOrdering();
  void writeMaxLineNum();
};

}  // namespace combgraph
}  // namespace octi

#endif  // OCTI_COMBGRAPH_GRAPH_H_
