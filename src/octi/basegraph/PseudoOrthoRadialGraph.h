// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_PSEUDOORTHORADIALGRAPH_H_
#define OCTI_BASEGRAPH_PSEUDOORTHORADIALGRAPH_H_

#include "octi/basegraph/GridGraph.h"

namespace octi {
namespace basegraph {

class PseudoOrthoRadialGraph : public GridGraph {
 public:
  using GridGraph::neigh;
  PseudoOrthoRadialGraph(const util::geo::DBox& bbox, double cellSize,
                         double spacer, const Penalties& pens)
      : GridGraph(util::geo::extendBox(
                      bbox, util::geo::getBoundingBox(util::geo::rotate(
                                convexHull(bbox), 90, centroid(bbox)))),
                  cellSize, spacer, pens) {
    _numBeams = 8;
    _heurHopCost = _c.p_45 - _c.p_135;
  }

  virtual void init();
  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
  getHeur(const std::set<GridNode*>& to) const;
  virtual double heurCost(int64_t xa, int64_t ya, int64_t xb, int64_t yb) const;

  virtual PolyLine<double> geomFromPath(
      const std::vector<std::pair<size_t, size_t>>& res) const;
  virtual double ndMovePen(const CombNode* cbNd, const GridNode* grNd) const;

 protected:
  virtual void writeInitialCosts();
  virtual GridNode* writeNd(size_t x, size_t y);
  virtual GridNode* neigh(size_t cx, size_t cy, size_t i) const;
  virtual GridNode* getNode(size_t x, size_t y) const;
  virtual void getSettledAdjEdgs(GridNode* n, CombNode* origNd,
                                 CombEdge* outgoing[8]);

 private:
  virtual int multi(size_t y) const;

  size_t _numBeams;
};

struct PseudoOrthoRadialGraphHeur
    : public util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float> {
  PseudoOrthoRadialGraphHeur(const basegraph::GridGraph* g,
                             const std::set<GridNode*>& to)
      : g(g), to(0) {
    cheapestSink = std::numeric_limits<float>::infinity();

    for (auto n : to) {
      assert(n->pl().getParent() == n);
      size_t i = 0;
      for (; i < g->maxDeg(); i++) {
        if (!n->pl().getPort(i)) continue;
        float sinkCost = g->getEdg(n->pl().getPort(i), n)->pl().cost();
        if (sinkCost < cheapestSink) cheapestSink = sinkCost;
        auto neigh = g->neigh(n, i);
        if (neigh && to.find(neigh) == to.end()) {
          hull.push_back(n->pl().getX());
          hull.push_back(n->pl().getY());
          break;
        }
      }
      for (size_t j = i; j < g->maxDeg(); j++) {
        if (!n->pl().getPort(j)) continue;
        float sinkCost = g->getEdg(n->pl().getPort(j), n)->pl().cost();
        if (sinkCost < cheapestSink) cheapestSink = sinkCost;
      }
    }
  }

  float operator()(const GridNode* from, const std::set<GridNode*>& to) const {
    if (to.count(from->pl().getParent())) return 0;

    float ret = std::numeric_limits<float>::infinity();

    for (size_t i = 0; i < hull.size(); i += 2) {
      float tmp = g->heurCost(from->pl().getParent()->pl().getX(),
                              from->pl().getParent()->pl().getY(), hull[i],
                              hull[i + 1]);
      if (tmp < ret) ret = tmp;
    }

    return ret + cheapestSink;
  }

  const octi::basegraph::BaseGraph* g;
  GridNode* to;
  std::vector<size_t> hull;
  float cheapestSink;
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_PSEUDOORTHORADIALGRAPH_H_
