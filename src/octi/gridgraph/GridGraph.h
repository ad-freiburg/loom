// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_GRIDGRAPH_H_
#define OCTI_GRIDGRAPH_GRIDGRAPH_H_

#include <queue>
#include <set>
#include <unordered_map>
#include "octi/combgraph/CombGraph.h"
#include "octi/gridgraph/GridEdgePL.h"
#include "octi/gridgraph/GridNodePL.h"
#include "octi/gridgraph/NodeCost.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/DirGraph.h"

#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"

using util::graph::DirGraph;
using util::graph::Node;
using util::geo::Grid;
using util::geo::Point;

using octi::combgraph::CombNode;
using octi::combgraph::CombEdge;

namespace octi {
namespace gridgraph {

typedef util::graph::Node<GridNodePL, GridEdgePL> GridNode;
typedef util::graph::Edge<GridNodePL, GridEdgePL> GridEdge;

typedef std::vector<double> GeoPens;
typedef std::map<const CombEdge*, GeoPens> GeoPensMap;

struct Candidate {
  Candidate(GridNode* n, double d) : n(n), d(d){};

  bool operator<(const Candidate& c) const { return d > c.d; }

  GridNode* n;
  double d;
};

struct Penalties {
  double p_0, p_45, p_90, p_135;
  double verticalPen, horizontalPen, diagonalPen;
  double densityPen;
};

class GridGraph : public DirGraph<GridNodePL, GridEdgePL> {
 public:
  GridGraph(const util::geo::DBox& bbox, double cellSize, double spacer,
            const Penalties& pens);

  double getCellSize() const;

  GridNode* getNode(size_t x, size_t y) const;

  const Grid<GridNode*, Point, double>& getGrid() const;

  NodeCost nodeBendPenalty(GridNode* n, CombEdge* e);
  NodeCost topoBlockPenalty(GridNode* n, CombNode* origNode, CombEdge* e);
  NodeCost spacingPenalty(GridNode* n, CombNode* origNode, CombEdge* e);

  double heurCost(int64_t xa, int64_t ya, int64_t xb, int64_t yb) const;

  std::priority_queue<Candidate> getGridNdCands(const util::geo::DPoint& p,
                                                double maxD) const;

  void addCostVector(GridNode* n, const NodeCost& addC);

  void openNodeSinkTo(GridNode* n, double cost);
  void closeNodeSinkTo(GridNode* n);
  void openNodeSinkFr(GridNode* n, double cost);
  void closeNodeSinkFr(GridNode* n);
  void openNodeTurns(GridNode* n);
  void closeNodeTurns(GridNode* n);

  GridNode* getNeighbor(size_t cx, size_t cy, size_t i) const;

  std::set<GridNode*> getGrNdCands(CombNode* n, double maxDis);

  void settleNd(GridNode* n, CombNode* cn);
  void settleEdg(GridNode* a, GridNode* b, CombEdge* e);

  const Penalties& getPenalties() const;

  void unSettleEdg(GridNode* a, GridNode* b);
  void unSettleNd(CombNode* a);

  bool isSettled(const CombNode* cn);

  GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  void reset();

  GridNode* getSettled(const CombNode* cnd) const;

  double ndMovePen(const CombNode* cbNd, const GridNode* grNd) const;

  GridNode* getGrNdById(size_t id) const;
  const GridEdge* getGrEdgById(std::pair<size_t, size_t> id) const;
  void addResEdg(GridEdge* ge, CombEdge* cg);
  CombEdge* getResEdg(GridEdge* ge);

  void writeGeoCoursePens(const CombEdge* ce, GeoPensMap* target, double pen);

 private:
  util::geo::DBox _bbox;
  Penalties _c;

  Grid<GridNode*, Point, double> _grid;
  double _cellSize, _spacer;
  std::unordered_map<const CombNode*, GridNode*> _settled;

  double _heurECost, _heurHopCost;

  // encoding portable IDs for each node
  std::vector<GridNode*> _nds;

  // edge id counter
  size_t _edgeCount;

  std::unordered_map<GridEdge*, CombEdge*> _resEdgs;

  void writeInitialCosts();

  GridNode* writeNd(size_t x, size_t y);

  void getSettledAdjEdgs(GridNode* n, CombEdge* outgoing[8]);
};
}
}

#endif  // OCTI_GRIDGRAPH_GRIDGRAPH_H_
