// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_GRIDGRAPH_H_
#define OCTI_BASEGRAPH_GRIDGRAPH_H_

#include <queue>
#include <set>
#include <unordered_map>
#include "octi/combgraph/CombGraph.h"
#include "octi/basegraph/GridEdgePL.h"
#include "octi/basegraph/GridNodePL.h"
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/BaseGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/DirGraph.h"

#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/CombNodePL.h"

using util::geo::Grid;
using util::geo::Point;
using util::graph::DirGraph;
using util::graph::Node;

using octi::combgraph::CombEdge;
using octi::combgraph::CombNode;

namespace octi {
namespace basegraph {

class GridGraph : public BaseGraph {
 public:
  GridGraph(const util::geo::DBox& bbox, double cellSize, double spacer,
            const Penalties& pens);

  virtual double getCellSize() const;

  virtual NodeCost nodeBendPen(GridNode* n, CombEdge* e);
  virtual NodeCost topoBlockPen(GridNode* n, CombNode* origNode, CombEdge* e);
  virtual NodeCost spacingPen(GridNode* n, CombNode* origNode, CombEdge* e);

  virtual double heurCost(int64_t xa, int64_t ya, int64_t xb, int64_t yb) const;

  virtual std::priority_queue<Candidate> getGridNdCands(const util::geo::DPoint& p,
                                                double maxD) const;

  virtual void addCostVec(GridNode* n, const NodeCost& addC);

  virtual void openSinkTo(GridNode* n, double cost);
  virtual void closeSinkTo(GridNode* n);
  virtual void openSinkFr(GridNode* n, double cost);
  virtual void closeSinkFr(GridNode* n);
  virtual void openTurns(GridNode* n);
  virtual void closeTurns(GridNode* n);

  virtual GridNode* getNeighbor(const GridNode* n, size_t i) const;
  virtual size_t getNumNeighbors() const;

  virtual std::set<GridNode*> getGrNdCands(CombNode* n, double maxDis);

  virtual void settleNd(GridNode* n, CombNode* cn);
  virtual void settleEdg(GridNode* a, GridNode* b, CombEdge* e);

  virtual const Penalties& getPens() const;

  virtual void unSettleEdg(GridNode* a, GridNode* b);
  virtual void unSettleNd(CombNode* a);

  virtual bool isSettled(const CombNode* cn);

  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual void reset();

  virtual GridNode* getSettled(const CombNode* cnd) const;

  virtual double ndMovePen(const CombNode* cbNd, const GridNode* grNd) const;

  virtual GridNode* getGrNdById(size_t id) const;
  virtual const GridEdge* getGrEdgById(std::pair<size_t, size_t> id) const;
  virtual void addResEdg(GridEdge* ge, CombEdge* cg);
  virtual CombEdge* getResEdg(GridEdge* ge);

  virtual CrossEdgPairs getCrossEdgPairs() const;

  virtual void writeGeoCoursePens(const CombEdge* ce, GeoPensMap* target, double pen);

  virtual void addObstacle(const util::geo::Polygon<double>& obst);

 protected:
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

  std::vector<util::geo::Polygon<double>> _obstacles;

  std::unordered_map<GridEdge*, CombEdge*> _resEdgs;

  const Grid<GridNode*, Point, double>& getGrid() const;

  void writeInitialCosts();
  void writeObstacleCost(const util::geo::Polygon<double>& obst);
  void reWriteObstCosts();
  void reWriteGeoPens();

  GridNode* getNode(size_t x, size_t y) const;

  GridNode* writeNd(size_t x, size_t y);

  GridNode* getNeighbor(size_t cx, size_t cy, size_t i) const;

  void getSettledAdjEdgs(GridNode* n, CombEdge* outgoing[8]);
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_GRIDGRAPH_H_
