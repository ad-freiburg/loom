// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_GRIDGRAPH_H_
#define OCTI_GRIDGRAPH_GRIDGRAPH_H_

#include <queue>
#include <unordered_set>
#include <unordered_map>
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/Graph.h"
#include "octi/gridgraph/NodePL.h"
#include "octi/gridgraph/EdgePL.h"

#include "octi/graph/CombNodePL.h"
#include "octi/graph/CombEdgePL.h"

using util::graph::Graph;
using util::graph::Node;
using util::geo::Grid;
using util::geo::Point;

namespace octi {
namespace gridgraph {

typedef util::graph::Node<gridgraph::NodePL, gridgraph::EdgePL> GridNode;
typedef util::graph::Edge<gridgraph::NodePL, gridgraph::EdgePL> GridEdge;

typedef util::graph::Graph<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombGraph;
typedef util::graph::Node<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombNode;
typedef util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>
    CombEdge;

struct Candidate {
  Candidate(Node<NodePL, EdgePL>* n, double d) : n(n), d(d) {};

  bool operator<(const Candidate& c) const { return d > c.d; }

  Node<NodePL, EdgePL>* n;
  double d;
};

struct Penalties {
  double p_0, p_45, p_90, p_135;
  double verticalPen, horizontalPen, diagonalPen;
};

class GridGraph : public Graph<NodePL, EdgePL> {
 public:
  GridGraph(const util::geo::Box& bbox, double cellSize, const Penalties& pens);

  Node<NodePL, EdgePL>* getNode(size_t x, size_t y) const;

  const Grid<Node<NodePL, EdgePL>*, Point>& getGrid() const;

  void spacingPenalty(GridNode* n, CombNode* origNode, CombEdge* e, double* ret);
  void topoBlockPenalty(GridNode* n, CombNode* origNode, CombEdge* e, double* ret);
  void outDegDeviationPenalty(GridNode* n, CombNode* origNode, CombEdge* e, double* addC);
  void balanceEdge(GridNode* a, GridNode* b);

  double heurCost(int64_t xa, int64_t ya, int64_t xb, int64_t yb) const;

  std::priority_queue<Candidate> getNearestCandidatesFor(const util::geo::Point& p, double maxD) const;

  void addCostVector(GridNode* n, double addC[8], double* invCost);
  void removeCostVector(GridNode* n, double addC[8]);
  std::pair<size_t, size_t> getNodeCoords(GridNode* n) const;

  void openNodeSink(GridNode* n, double cost);
  void closeNodeSink(GridNode* n);
  void openNode(GridNode* n);
  void closeNode(GridNode* n);

  Node<NodePL, EdgePL>* getNeighbor(size_t cx, size_t cy, size_t i) const;

  GridNode* getGridNodeFrom(CombNode* n, double maxDis);
  std::unordered_set<GridNode*> getGridNodesTo(CombNode* n, double maxDis);

  void settleGridNode(GridNode* n, CombNode* cn);
  bool isSettled(CombNode* cn);
 private:
  util::geo::Box _bbox;
  Penalties _c;

  Grid<Node<NodePL, EdgePL>*, Point> _grid;
  std::unordered_map<CombNode*, GridNode*> _settled;

  std::set<util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>*> getResEdges(Node<NodePL, EdgePL>* n) const;

  void writeInitialCosts();

  GridEdge* getNEdge(GridNode* a, GridNode* b);
  GridEdge* getOtherEdge(GridEdge* e);
  void getSettledOutgoingEdges(GridNode* n, CombEdge* outgoing[8]);
};

}
}

#endif  // OCTI_GRIDGRAPH_GRIDGRAPH_H_
