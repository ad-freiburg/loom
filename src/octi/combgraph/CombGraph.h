// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_GRAPH_H_
#define OCTI_COMBGRAPH_GRAPH_H_

#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"
#include "octi/transitgraph/TransitGraph.h"
#include "util/graph/UndirGraph.h"

namespace octi {
namespace combgraph {

typedef util::graph::Node<CombNodePL, CombEdgePL> CombNode;
typedef util::graph::Edge<CombNodePL, CombEdgePL> CombEdge;

using octi::transitgraph::TransitGraph;
using octi::transitgraph::EdgeOrdering;

class CombGraph : public util::graph::UndirGraph<CombNodePL, CombEdgePL> {
 public:
  CombGraph(const TransitGraph* g);
  void getTransitGraph(TransitGraph* target) const;

  EdgeOrdering getEdgeOrderingForNode(CombNode* n) const;
  EdgeOrdering getEdgeOrderingForNode(
      CombNode* n, bool useOrigNextNode,
      const std::map<CombNode*, DPoint>& newPos) const;
  size_t changesTopology(CombNode* nOrig, DPoint p,
                         const std::map<CombNode*, DPoint>& newPos) const;

 private:
  void build(const TransitGraph* source);
  void combineDeg2();
  void writeEdgeOrdering();
};

}  // combgraph
}  // octi

#endif  // OCTI_COMBGRAPH_GRAPH_H_
