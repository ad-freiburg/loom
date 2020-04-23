#ifndef OCTI_COMBGRAPH_DRAWING_H_
#define OCTI_COMBGRAPH_DRAWING_H_

#include <map>
#include "octi/combgraph/CombGraph.h"
#include "octi/basegraph/BaseGraph.h"
#include "util/graph/Dijkstra.h"

namespace octi {
namespace combgraph {

using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::CombEdge;
using util::graph::Dijkstra;
using octi::basegraph::BaseGraph;
using octi::basegraph::GridNode;
using octi::basegraph::GridEdge;
using octi::basegraph::GridNodePL;
using octi::basegraph::GridEdgePL;

typedef Dijkstra::EList<GridNodePL, GridEdgePL> GrEdgList;

struct Costs {
  double bend;
  double move;
  double hop;
  double dense;
};

class Drawing {
 public:
  Drawing(const BaseGraph* gg)
      : _c(std::numeric_limits<double>::infinity()), _gg(gg){};
  Drawing() : _c(std::numeric_limits<double>::infinity()), _gg(0){};

  double score() const;
  Costs fullScore() const;
  void crumble();

  void draw(CombEdge* ce, const GrEdgList& ge, bool rev);
  void erase(CombEdge* ce);
  void erase(CombNode* ce);

  void getLineGraph(LineGraph* target) const;

  const GridNode* getGrNd(const CombNode* cn);

  bool drawn(const CombEdge* ce) const;

  void eraseFromGrid(const CombEdge* ce, BaseGraph* gg);
  void eraseFromGrid(const CombNode* ce, BaseGraph* gg);
  void applyToGrid(const CombEdge* ce, BaseGraph* gg);
  void applyToGrid(const CombNode* ce, BaseGraph* gg);

  void eraseFromGrid(BaseGraph* gg);
  void applyToGrid(BaseGraph* gg);

  double getEdgCost(const CombEdge* e) const;
  double getNdBndCost(const CombNode* e) const;
  double getNdReachCost(const CombNode* e) const;

  void setBaseGraph(const BaseGraph* gg);

 private:
  std::map<const CombNode*, size_t> _nds;
  std::map<const CombEdge*, std::vector<std::pair<size_t, size_t>>> _edgs;

  std::map<const CombNode*, double> _ndReachCosts;
  std::map<const CombNode*, double> _ndBndCosts;
  std::map<const CombEdge*, double> _edgCosts;
  std::map<const CombEdge*, double> _springCosts;
  double _c;
  const BaseGraph* _gg;

  double recalcBends(const CombNode* nd);
};
}
}

#endif  // OCTI_COMBGRAPH_DRAWING_H_
