// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_BASEGRAPH_H_
#define OCTI_BASEGRAPH_BASEGRAPH_H_

#include <queue>
#include <set>
#include <unordered_map>
#include "octi/basegraph/GridEdgePL.h"
#include "octi/basegraph/GridNodePL.h"
#include "octi/basegraph/NodeCost.h"
#include "octi/combgraph/CombGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/DirGraph.h"
#include "util/graph/Dijkstra.h"

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

typedef util::graph::Node<GridNodePL, GridEdgePL> GridNode;
typedef util::graph::Edge<GridNodePL, GridEdgePL> GridEdge;

typedef std::pair<const GridEdge*, const GridEdge*> EdgPair;
typedef std::vector<std::pair<EdgPair, EdgPair>> CrossEdgPairs;

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

class BaseGraph : public DirGraph<GridNodePL, GridEdgePL> {
 public:
  BaseGraph(){};

  virtual void init() = 0;
  virtual double getCellSize() const = 0;

  virtual NodeCost nodeBendPen(GridNode* n, CombEdge* e) = 0;
  virtual NodeCost topoBlockPen(GridNode* n, CombNode* origNode,
                                CombEdge* e) = 0;
  virtual NodeCost spacingPen(GridNode* n, CombNode* origNode, CombEdge* e) = 0;

  virtual double heurCost(int64_t xa, int64_t ya, int64_t xb,
                          int64_t yb) const = 0;

  virtual const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
  getHeur(const std::set<GridNode*>& to) const = 0;

  virtual std::priority_queue<Candidate> getGridNdCands(
      const util::geo::DPoint& p, double maxD) const = 0;

  virtual void addCostVec(GridNode* n, const NodeCost& addC) = 0;

  virtual void openSinkTo(GridNode* n, double cost) = 0;
  virtual void closeSinkTo(GridNode* n) = 0;
  virtual void openSinkFr(GridNode* n, double cost) = 0;
  virtual void closeSinkFr(GridNode* n) = 0;
  virtual void openTurns(GridNode* n) = 0;
  virtual void closeTurns(GridNode* n) = 0;

  virtual GridNode* neigh(const GridNode* n, size_t i) const = 0;

  virtual size_t maxDeg() const = 0;

  virtual std::set<GridNode*> getGrNdCands(CombNode* n, double maxDis) = 0;

  virtual void settleNd(GridNode* n, CombNode* cn) = 0;
  virtual void settleEdg(GridNode* a, GridNode* b, CombEdge* e) = 0;

  virtual const Penalties& getPens() const = 0;

  virtual void unSettleEdg(GridNode* a, GridNode* b) = 0;
  virtual void unSettleNd(CombNode* a) = 0;

  virtual bool isSettled(const CombNode* cn) = 0;

  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const = 0;
  virtual void reset() = 0;

  virtual GridNode* getSettled(const CombNode* cnd) const = 0;

  virtual double ndMovePen(const CombNode* cbNd,
                           const GridNode* grNd) const = 0;

  virtual double getBendPen(size_t origI, size_t targetI) const = 0;
  virtual GridNode* getGrNdById(size_t id) const = 0;
  virtual const GridEdge* getGrEdgById(std::pair<size_t, size_t> id) const = 0;
  virtual void addResEdg(GridEdge* ge, CombEdge* cg) = 0;
  virtual CombEdge* getResEdg(GridEdge* ge) = 0;

  virtual void writeGeoCoursePens(const CombEdge* ce, GeoPensMap* target,
                                  double pen) = 0;

  virtual CrossEdgPairs getCrossEdgPairs() const = 0;

  virtual void addObstacle(const util::geo::Polygon<double>& obst) = 0;
  virtual PolyLine<double> geomFromPath(const std::vector<std::pair<size_t, size_t>>& res) const = 0;
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_BASEGRAPH_H_
