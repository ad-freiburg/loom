#ifndef OCTI_COMBGRAPH_DRAWING_H_
#define OCTI_COMBGRAPH_DRAWING_H_

#include <map>
#include "octi/basegraph/BaseGraph.h"
#include "octi/combgraph/CombGraph.h"
#include "util/graph/Dijkstra.h"

namespace octi {
namespace combgraph {

using octi::basegraph::BaseGraph;
using octi::basegraph::GridEdge;
using octi::basegraph::GridEdgePL;
using octi::basegraph::GridNode;
using octi::basegraph::GridNodePL;
using octi::combgraph::CombEdge;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using util::graph::Dijkstra;

typedef Dijkstra::EList<GridNodePL, GridEdgePL> GrEdgList;
typedef std::vector<std::pair<size_t, size_t>>  GrPath;

struct Score {
  Score(double bend, double move, double hop, double dense, double full, uint64_t violations) : bend(bend), move(move), hop(hop), dense(dense), full(full), violations(violations) {}
  Score() : bend(0), move(0), hop(0), dense(0), full(0), violations(0) {}
  double bend;
  double move;
  double hop;
  double dense;
  double full;
  uint64_t violations;
};

struct NodeOnSeg {
  shared::linegraph::LineNode* n;
  double progr;
};

inline bool operator<(const NodeOnSeg& lh, const NodeOnSeg& rh) {
  return lh.progr < rh.progr;
}

struct LineEdgeOnSeg {
  shared::linegraph::LineEdge* n;
  shared::linegraph::LineNode* start;
  shared::linegraph::LineNode* end;
};

struct Segment {
  GridNode* start;
  GridNode* end;
  std::vector<std::pair<size_t, size_t>> path;
  PolyLine<double> geom;
  std::vector<NodeOnSeg> nodes;
  std::set<CombEdge*> combEdges;
};

class Drawing {
 public:
  Drawing(const BaseGraph* gg)
      : _c(std::numeric_limits<double>::infinity()), _gg(gg), _violations(0){};
  Drawing() : _c(std::numeric_limits<double>::infinity()), _gg(0), _violations(0){};

  double score() const;
  double rawScore() const;
  uint64_t violations() const;
  Score fullScore() const;
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

  const std::map<const CombEdge*, GrPath>& getEdgPaths() const;

 private:
  std::map<const CombNode*, size_t> _nds;
  std::map<const CombEdge*, GrPath> _edgs;

  std::map<const CombNode*, double> _ndReachCosts;
  std::map<const CombNode*, double> _ndBndCosts;
  std::map<const CombEdge*, double> _edgCosts;
  std::map<const CombEdge*, int> _vios;
  std::map<const CombEdge*, double> _springCosts;
  double _c;
  const BaseGraph* _gg;

  size_t _violations;

  double recalcBends(const CombNode* nd);
};
}  // namespace combgraph
}  // namespace octi

#endif  // OCTI_COMBGRAPH_DRAWING_H_
