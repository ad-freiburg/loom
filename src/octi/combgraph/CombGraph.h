// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_GRAPH_H_
#define OCTI_COMBGRAPH_GRAPH_H_

#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"
#include "util/graph/UndirGraph.h"
#include "octi/transitgraph/TransitGraph.h"

namespace octi {
namespace combgraph {

typedef util::graph::Node<CombNodePL, CombEdgePL> CombNode;
typedef util::graph::Edge<CombNodePL, CombEdgePL> CombEdge;

using octi::transitgraph::TransitGraph;

class CombGraph : public util::graph::UndirGraph<CombNodePL, CombEdgePL> {
 public:
   CombGraph(const TransitGraph* g);
   void getTransitGraph(TransitGraph* target) const;
 private:
  void build(const TransitGraph* source);
  void combineDeg2();
};

}  // combgraph
}  // octi

#endif  // OCTI_COMBGRAPH_GRAPH_H_
