#ifndef OCTI_COMBGRAPH_DRAWING_H_
#define OCTI_COMBGRAPH_DRAWING_H_

#include <map>
#include "octi/combgraph/CombGraph.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/graph/Dijkstra.h"

namespace octi {
namespace combgraph {

using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::CombEdge;
using util::graph::Dijkstra;
using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;
using octi::gridgraph::GridNodePL;
using octi::gridgraph::GridEdgePL;

typedef Dijkstra::EList<GridNodePL, GridEdgePL> GrEdgList;

class Drawing {
 public:
  Drawing() : _c(0){};
  Drawing(double cost) : _c(cost){};

  double score() const;

  void draw(CombEdge* ce, const GrEdgList& ge);

  void getTransitGraph(TransitGraph* target) const;

 private:
  std::map<const CombNode*, GridNode*> _nds;
  std::map<const CombEdge*, std::vector<GridEdge*>> _edgs;
  double _c;

PolyLine<double> buildPolylineFromRes(
    const std::vector<GridEdge*>& res) const;
};
}
}

#endif  // OCTI_COMBGRAPH_DRAWING_H_
