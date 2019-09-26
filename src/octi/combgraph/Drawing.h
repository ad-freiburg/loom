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
  Drawing(const GridGraph* gg) : _c(0), _gg(gg) {};

  double score() const;

  void draw(CombEdge* ce, const GrEdgList& ge, bool rev);

  void getTransitGraph(TransitGraph* target) const;

  const GridNode* getGrNd(const CombNode* cn);
  const std::vector<const GridEdge*>& getGrEdgs(const CombEdge* ce);

 private:
  std::map<const CombNode*, const GridNode*> _nds;
  std::map<const CombEdge*, std::vector<const GridEdge*>> _edgs;
  double _c;
  const GridGraph* _gg;

PolyLine<double> buildPolylineFromRes(
    const std::vector<const GridEdge*>& res) const;
};
}
}

#endif  // OCTI_COMBGRAPH_DRAWING_H_
