// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_OCTILINEARIZER_H_
#define OCTI_OCTILINEARIZER_H_

#include <unordered_set>
#include <vector>
#include "octi/basegraph/BaseGraph.h"
#include "octi/basegraph/GridGraph.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
#include "shared/linegraph/LineGraph.h"
#include "util/graph/BiDijkstra.h"
#include "util/graph/Dijkstra.h"

namespace octi {

using octi::basegraph::GeoPens;
using octi::basegraph::GeoPensMap;
using octi::basegraph::GridEdge;
using octi::basegraph::GridEdgePL;
using octi::basegraph::GridGraph;
using octi::basegraph::GridNode;
using octi::basegraph::GridNodePL;
using octi::basegraph::NodeCost;
using octi::basegraph::Penalties;

using shared::linegraph::LineEdge;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;

using octi::combgraph::CombEdge;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;
using octi::combgraph::Drawing;
using octi::combgraph::EdgeOrdering;

using octi::combgraph::Score;

using util::graph::Dijkstra;

typedef util::graph::EList<GridNodePL, GridEdgePL> GrEdgList;
typedef util::graph::NList<GridNodePL, GridEdgePL> GrNdList;
typedef std::pair<std::set<GridNode*>, std::set<GridNode*>> RtPair;
typedef std::map<CombNode*, const GridNode*> SettledPos;

enum Undrawable { DRAWN = 0, NO_PATH = 1, NO_CANDS = 2 };

// exception thrown when no planar embedding could be found
struct NoEmbeddingFoundExc : public std::exception {
  const char* what() const throw() {
    return "Could not find planar embedding for input graph.";
  }
};

// comparator for nodes, based on degree
struct NodeCmp {
  bool operator()(CombNode* a, CombNode* b) {
    return a->pl().getRouteNumber() < b->pl().getRouteNumber();
  }
};

typedef std::priority_queue<CombNode*, std::vector<CombNode*>, NodeCmp> NodePQ;

struct GraphMeasures {
  double maxNodeDist;
  double minNodeDist;
  double avgNodeDist;
  double maxEdgeLength;
  double minEdgeLength;
  double avgEdgeLength;
  size_t maxDeg;
};

struct GridCost
    : public util::graph::Dijkstra::CostFunc<GridNodePL, GridEdgePL, float> {
  GridCost(float inf) : _inf(inf) {}
  virtual float operator()(const GridNode* from, const GridEdge* e,
                           const GridNode* to) const {
    UNUSED(from);
    UNUSED(to);
    return e->pl().cost();
  }

  float _inf;

  virtual float inf() const { return _inf; }
};

struct GridCostGeoPen
    : public Dijkstra::CostFunc<GridNodePL, GridEdgePL, float> {
  GridCostGeoPen(float inf, const GeoPens* geoPens)
      : _inf(inf), _geoPens(geoPens) {}
  virtual float operator()(const GridNode* from, const GridEdge* e,
                           const GridNode* to) const {
    UNUSED(from);
    UNUSED(to);
    return e->pl().cost() + (*_geoPens)[e->pl().getId()];
  }

  float _inf;
  const GeoPens* _geoPens;

  virtual float inf() const { return _inf; }
};

class Octilinearizer {
 public:
  Octilinearizer(basegraph::BaseGraphType baseGraphType)
      : _baseGraphType(baseGraphType) {}

  Score draw(const CombGraph& cg, const util::geo::DBox& box, LineGraph* out,
             basegraph::BaseGraph** gg, const Penalties& pens, double gridSize,
             double borderRad, double maxGrDist, bool restrLocSearch,
             double enfGeoCourse,
             const std::vector<util::geo::Polygon<double>>& obstacles,
             size_t abortAfter);

  Score drawILP(const CombGraph& cg, const util::geo::DBox& box, LineGraph* out,
                basegraph::BaseGraph** gg, const Penalties& pens,
                double gridSize, double borderRad, double maxGrDist,
                bool noSolve, double enfGeoPens, int timeLim,
                const std::string& solverStr, const std::string& path);

 private:
  basegraph::BaseGraphType _baseGraphType;

  basegraph::BaseGraph* newBaseGraph(const util::geo::DBox& bbox,
                                     const CombGraph& cg, double cellSize,
                                     double spacer,
                                     const Penalties& pens) const;

  void writeNdCosts(GridNode* n, CombNode* origNode, CombEdge* e,
                    basegraph::BaseGraph* g);

  void settleRes(GridNode* startGridNd, GridNode* toGridNd,
                 basegraph::BaseGraph* gg, CombNode* from, CombNode* to,
                 const GrEdgList& res, CombEdge* e, size_t rndrOrder);

  static const CombNode* getCenterNd(const CombGraph* cg);

  std::vector<CombEdge*> getOrdering(const CombGraph& cg, bool randr) const;

  Undrawable draw(const std::vector<CombEdge*>& order, basegraph::BaseGraph* gg,
                  Drawing* drawing, double cutoff, double maxGrDist,
                  const GeoPensMap* geoPensMap, size_t abortAfter);
  Undrawable draw(const std::vector<CombEdge*>& order,
                  const SettledPos& settled, basegraph::BaseGraph* gg,
                  Drawing* drawing, double cutoff, double maxGrDist,
                  const GeoPensMap* geoPensMap, size_t abortAfter);

  SettledPos neigh(const SettledPos& pos, const std::vector<CombNode*>&,
                   size_t i) const;

  RtPair getRtPair(CombNode* frCmbNd, CombNode* toCmbNd,
                   const SettledPos& settled, basegraph::BaseGraph* gg,
                   double maxGrDist);

  std::set<GridNode*> getCands(CombNode* cmBnd, const SettledPos& settled,
                               basegraph::BaseGraph* gg, double maxDis);
};

}  // namespace octi

#endif  // OCTI_OCTILINEARIZER_H_
