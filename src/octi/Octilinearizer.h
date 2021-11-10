// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_OCTILINEARIZER_H_
#define OCTI_OCTILINEARIZER_H_

#include <unordered_set>
#include <vector>
#include "ilp/ILPGridOptimizer.h"
#include "octi/basegraph/BaseGraph.h"
#include "octi/basegraph/GridGraph.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
#include "octi/config/OctiConfig.h"
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
struct NodeCmpDeg {
  bool operator()(const CombNode* a, const CombNode* b) {
    // smallest first, as the PQ returns the biggest
    return a->getDeg() < b->getDeg();
  }
};

// comparator for nodes, based on line degree
struct NodeCmpLdeg {
  bool operator()(const CombNode* a, const CombNode* b) {
    // smallest first, as the PQ returns the biggest
    return a->pl().getLDeg() < b->pl().getLDeg();
  }
};

struct EdgeCmpNumLines {
  bool operator()(const CombEdge* a, const CombEdge* b) {
    return a->pl().getNumLines() > b->pl().getNumLines();
  }
};

struct EdgeCmpLength {
  bool operator()(const CombEdge* a, const CombEdge* b) {
    double da = util::geo::dist(*a->getFrom()->pl().getGeom(),
                                *a->getTo()->pl().getGeom());
    double db = util::geo::dist(*b->getFrom()->pl().getGeom(),
                                *b->getTo()->pl().getGeom());

    // smallest first, as we are sorting a vector with this
    return da < db;
  }
};

struct EdgeCmpLdeg {
  bool operator()(const CombEdge* a, const CombEdge* b) {
    std::pair<size_t, size_t> oa = {
        std::max(a->getFrom()->pl().getLDeg(), a->getTo()->pl().getLDeg()),
        std::min(a->getFrom()->pl().getLDeg(), a->getTo()->pl().getLDeg())};
    std::pair<size_t, size_t> ob = {
        std::max(b->getFrom()->pl().getLDeg(), b->getTo()->pl().getLDeg()),
        std::min(b->getFrom()->pl().getLDeg(), b->getTo()->pl().getLDeg())};

    // biggest first, as we are sorting a vector with this
    return oa > ob;
  }
};

struct EdgeCmpDeg {
  bool operator()(const CombEdge* a, const CombEdge* b) {
    std::pair<size_t, size_t> oa = {
        std::max(a->getFrom()->getDeg(), a->getTo()->getDeg()),
        std::min(a->getFrom()->getDeg(), a->getTo()->getDeg())};
    std::pair<size_t, size_t> ob = {
        std::max(b->getFrom()->getDeg(), b->getTo()->getDeg()),
        std::min(b->getFrom()->getDeg(), b->getTo()->getDeg())};

    // biggest first, as we are sorting a vector with this
    return oa > ob;
  }
};

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
             basegraph::BaseGraph** gg, Drawing* d, const Penalties& pens,
             double gridSize, double borderRad, double maxGrDist,
             config::OrderMethod orderMethod, bool restrLocSearch,
             double enfGeoCourse, size_t hananIters,
             const std::vector<util::geo::Polygon<double>>& obstacles,
             size_t locsearchIters, size_t abortAfter);

  Score drawILP(const CombGraph& cg, const util::geo::DBox& box, LineGraph* out,
                basegraph::BaseGraph** gg, Drawing* d, const Penalties& pens,
                double gridSize, double borderRad, double maxGrDist,
                config::OrderMethod orderMethod, bool noSolve,
                double enfGeoPens, size_t hananIters, int timeLim,
                const std::string& cacheDir, octi::ilp::ILPStats* stats,
                const std::string& solverStr, const std::string& path);

  size_t maxNodeDeg() const;

 private:
  basegraph::BaseGraphType _baseGraphType;

  basegraph::BaseGraph* newBaseGraph(const util::geo::DBox& bbox,
                                     const CombGraph& cg, double cellSize,
                                     double spacer, size_t hananIters,
                                     const Penalties& pens) const;

  util::geo::Polygon<double> hull(const CombGraph& cg) const;

  void writeNdCosts(GridNode* n, CombNode* origNode, CombEdge* e,
                    basegraph::BaseGraph* g);

  void settleRes(GridNode* startGridNd, GridNode* toGridNd,
                 basegraph::BaseGraph* gg, CombNode* from, CombNode* to,
                 const GrEdgList& res, CombEdge* e, size_t rndrOrder);

  static const CombNode* getCenterNd(const CombGraph* cg);

  std::vector<CombEdge*> getOrdering(const CombGraph& cg,
                                     octi::config::OrderMethod method) const;

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
                               basegraph::BaseGraph* gg, size_t maxGridDis);

  void statLine(Undrawable status, const std::string& msg,
                const Drawing& drawing, double ms,
                const std::string& mark) const;

  template <class ECmp, class NCmp>
  std::vector<CombEdge*> getGrowthOrder(const CombGraph& cg) const {
    ECmp eCmp;

    std::priority_queue<CombNode*, std::vector<CombNode*>, NCmp> globalPq,
        dangling;

    std::set<CombNode*> settled;
    std::vector<CombEdge*> retOrder;

    // global PQ in case graph is not connected
    for (auto n : cg.getNds()) globalPq.push(n);

    while (!globalPq.empty()) {
      auto n = globalPq.top();
      globalPq.pop();
      dangling.push(n);

      while (!dangling.empty()) {
        auto n = dangling.top();
        dangling.pop();

        if (settled.find(n) != settled.end()) continue;

        auto ordered = n->getAdjList();
        std::sort(ordered.begin(), ordered.end(), eCmp);

        for (auto ee : ordered) {
          if (settled.find(ee->getOtherNd(n)) != settled.end()) continue;
          dangling.push(ee->getOtherNd(n));
          retOrder.push_back(ee);
        }
        settled.insert(n);
      }
    }

    return retOrder;
  }
};

}  // namespace octi

#endif  // OCTI_OCTILINEARIZER_H_
