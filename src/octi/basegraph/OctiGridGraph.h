// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_OCTIGRIDGRAPH_H_
#define OCTI_BASEGRAPH_OCTIGRIDGRAPH_H_

#include "octi/basegraph/GridGraph.h"

namespace octi {
namespace basegraph {

class OctiGridGraph : public GridGraph {
 public:
  using GridGraph::getNeighbor;
  OctiGridGraph(const util::geo::DBox& bbox, double cellSize, double spacer,
                const Penalties& pens)
      : GridGraph(bbox, cellSize, spacer, pens) {

    _bendCosts[0] = _c.p_45 - _c.p_135;
    _bendCosts[3] = _c.p_45;
    _bendCosts[2] = _c.p_45 - _c.p_135 + _c.p_90;
    _bendCosts[1] = _bendCosts[0] + _bendCosts[3];

    // prepare the octigrid A* heuristic
    _heurXCost = _c.horizontalPen + _heurHopCost;
    _heurYCost = _c.verticalPen + _heurHopCost;
    _heurDiagCost = _c.diagonalPen + _heurHopCost;

    if (_heurDiagCost < _heurXCost) {
      // 1) two horizontal hops may be substituted by two diagonal hops
      // 2) a horizontal hop may be subsituted by a diagonal + a vertical hop
      // in both cases, the minimum x hops required is _heurDiagCost
      _heurXCost = _heurDiagCost;
    }

    if (_heurDiagCost < _heurYCost) {
      // 1) two vertical hops may be substituted by two diagonal hops
      // 2) a vertical hop may be subsituted by a diagonal + a horizontal hop
      // in both cases, the minimum y hops required is _heurDiagCost
      _heurYCost = _heurDiagCost;
    }

    if (_heurDiagCost > _heurXCost + _heurYCost) {
      // if the diagonal cost is greater than the horizontal plus the
      // vertical, we can never save something by going diagonal
      _heurDiagSave = 0;
    } else {
      // else, we save one x hop and one y hop at the cost of one
      // diagonal hop
      _heurDiagSave = _heurDiagCost - _heurXCost - _heurYCost;
    }
  }

  virtual void unSettleEdg(GridNode* a, GridNode* b);
  virtual void settleEdg(GridNode* a, GridNode* b, CombEdge* e);
  virtual CrossEdgPairs getCrossEdgPairs() const;
  virtual GridEdge* getNEdg(const GridNode* a, const GridNode* b) const;
  virtual size_t getNumNeighbors() const;

 protected:
  virtual void writeInitialCosts();
  virtual GridNode* writeNd(size_t x, size_t y);
  virtual GridNode* getNeighbor(size_t cx, size_t cy, size_t i) const;
  virtual GridNode* getNode(size_t x, size_t y) const;
  virtual double getBendPen(size_t origI, size_t targetI) const;
  virtual double heurCost(int64_t xa, int64_t ya, int64_t xb, int64_t yb) const;

 private:
  double _heurDiagSave;
  double _heurXCost;
  double _heurYCost;
  double _heurDiagCost;

  double _bendCosts[4];
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_OCTIGRIDGRAPH_H_
