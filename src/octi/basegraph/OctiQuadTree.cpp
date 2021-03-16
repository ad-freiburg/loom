// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/OctiQuadTree.h"
#include "util/Misc.h"
#include "util/geo/QuadTree.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::NodeCost;
using octi::basegraph::OctiQuadTree;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;
using util::geo::QuadNode;
using util::geo::QuadValue;
using util::geo::QuadTree;

// _____________________________________________________________________________
void OctiQuadTree::unSettleEdg(CombEdge* ce, GridNode* a, GridNode* b) {
  if (a == b) return;

  auto ge = getNEdg(a, b);
  auto gf = getNEdg(b, a);

  assert(ge);
  assert(gf);

  ge->pl().delResEdg();
  gf->pl().delResEdg();

  _resEdgs[ge].erase(ce);
  _resEdgs[gf].erase(ce);


  if (_resEdgs[ge].size() == 0) {
    if (!a->pl().isSettled()) openTurns(a);
    if (!b->pl().isSettled()) openTurns(b);
  }

  // unblock diagonal edges crossing this edge

  // for diagonal edges, x and y difference is equivalent
  int len = labs((int)a->pl().getX() - (int)b->pl().getX());

  size_t dir = getDir(a, b);

  GridNode* aa = 0;
  GridNode* bb = 0;

  if (dir == 1) {
    aa = getNode(a->pl().getX(), a->pl().getY() + len);
    bb = getNode(a->pl().getX() + len, a->pl().getY());
  } else if (dir == 3) {
    aa = getNode(a->pl().getX() + len, a->pl().getY());
    bb = getNode(a->pl().getX(), a->pl().getY() - len);
  } else if (dir == 5) {
    aa = getNode(a->pl().getX() - len, a->pl().getY());
    bb = getNode(a->pl().getX(), a->pl().getY() - len);
  } else if (dir == 7) {
    aa = getNode(a->pl().getX() - len, a->pl().getY());
    bb = getNode(a->pl().getX(), a->pl().getY() + len);
  }

  if (aa && bb) {
    auto e = getNEdg(aa, bb);
    auto f = getNEdg(bb, aa);
    if (e && f) {
      e->pl().unblock();
      f->pl().unblock();
    }
  }
}

// _____________________________________________________________________________
void OctiQuadTree::settleEdg(GridNode* a, GridNode* b, CombEdge* e,
                             size_t rndrOrd) {
  if (a == b) return;

  // this closes the grid edge
  auto ge = getNEdg(a, b);

  addResEdg(ge, e);

  ge->pl().setRndrOrder(rndrOrd);

  // this closes both nodes
  // a close means that all major edges reaching this node are closed
  closeTurns(a);
  closeTurns(b);

  // block diagonal edges crossing this edge
  // for diagonal edges, x and y difference is equivalent
  int len = labs((int)a->pl().getX() - (int)b->pl().getX());

  size_t dir = getDir(a, b);

  GridNode* aa = 0;
  GridNode* bb = 0;

  if (dir == 1) {
    aa = getNode(a->pl().getX(), a->pl().getY() + len);
    bb = getNode(a->pl().getX() + len, a->pl().getY());
  } else if (dir == 3) {
    aa = getNode(a->pl().getX() + len, a->pl().getY());
    bb = getNode(a->pl().getX(), a->pl().getY() - len);
  } else if (dir == 5) {
    aa = getNode(a->pl().getX() - len, a->pl().getY());
    bb = getNode(a->pl().getX(), a->pl().getY() - len);
  } else if (dir == 7) {
    aa = getNode(a->pl().getX() - len, a->pl().getY());
    bb = getNode(a->pl().getX(), a->pl().getY() + len);
  }

  if (aa && bb) {
    auto e = getNEdg(aa, bb);
    auto f = getNEdg(bb, aa);

    if (e && f) {
      e->pl().block();
      f->pl().block();
    }
  }
}

// _____________________________________________________________________________
CrossEdgPairs OctiQuadTree::getCrossEdgPairs() const {
  assert(false);
  CrossEdgPairs ret;
  for (const GridNode* n : getNds()) {
    if (!n->pl().isSink()) continue;

    auto eOr = getNEdg(n, neigh(n, 3));
    auto fOr = getNEdg(neigh(n, 3), n);

    if (!eOr || !fOr) continue;

    auto na = neigh(n, (3 + 7) % 8);
    auto nb = neigh(n, (3 + 1) % 8);

    if (!na || !nb) continue;

    auto e = getNEdg(na, nb);
    auto f = getNEdg(nb, na);
    ret.push_back({{eOr, fOr}, {e, f}});
  }

  return ret;
}

// _____________________________________________________________________________
void OctiQuadTree::init() {
  auto newBox = DBox();

  auto centroid = util::geo::centroid(_bbox);
  newBox = extendBox(_bbox, newBox);
  newBox = extendBox(rotate(convexHull(_bbox), 180, centroid), newBox);
  _bbox = newBox;

  size_t numCells =
      ceil(_bbox.getUpperRight().getX() - _bbox.getLowerLeft().getX()) /
      _cellSize;
  size_t n = pow(2, ceil(log2(numCells)));

  std::cerr << "num cells: " << n << std::endl;
  newBox = extendBox(DPoint(newBox.getLowerLeft().getX() + n * _cellSize,
                            newBox.getLowerLeft().getY() + n * _cellSize),
                     newBox);

  size_t maxDepth = log2(n);

  std::cerr << "Max depth: " << maxDepth << std::endl;

  _bbox = newBox;

struct SplitFunc : util::geo::SplitFunc<const CombNode*, double> {
  SplitFunc(double cellSize) : _cellSize(cellSize) {}
  virtual bool operator()(const QuadNode<double>& nd,
                       const QuadValue<const CombNode*, double>& newVal) const {
    UNUSED(newVal);
    double l = nd.bbox.getUpperRight().getX() - nd.bbox.getLowerLeft().getX();

    return nd.numEls > 0 || l > _cellSize;
  }
  double _cellSize;
} sFunc(_cellSize);

  QuadTree<const CombNode*, double> qt(maxDepth, sFunc, newBox);

  _grid = Grid<GridNode*, Point, double>(
      _cellSize, _cellSize, util::geo::pad(_bbox, _cellSize), false);

  // write nodes to quadtree
  for (auto cNd : _cg.getNds()) {
    qt.insert(cNd, *cNd->pl().getGeom());
  }

  struct SortByCellSize {
    SortByCellSize(const std::vector<QuadNode<double>>& qnds) : qnds(qnds) {}
    bool operator()(size_t aNid, size_t bNid) {
      const auto& a = qnds[aNid];
      const auto& b = qnds[bNid];
      return (a.bbox.getUpperRight().getX() - a.bbox.getLowerLeft().getX()) <
             (b.bbox.getUpperRight().getX() - b.bbox.getLowerLeft().getX());
    }
    const std::vector<QuadNode<double>>& qnds;
  } sortByCellSize(qt.getNds());

  std::vector<size_t> sortedQdNds;
  for (size_t i = 0; i < qt.getNds().size(); i++) {
    if (qt.getNd(i).numEls == -1) continue;
    sortedQdNds.push_back(i);
  }
  std::sort(sortedQdNds.begin(), sortedQdNds.end(), sortByCellSize);

  for (auto nid : sortedQdNds) {
    const auto& qNd = qt.getNd(nid);
    size_t xa = std::round(
        (qNd.bbox.getLowerLeft().getX() - _bbox.getLowerLeft().getX()) /
        _cellSize);
    size_t ya = std::round(
        (qNd.bbox.getLowerLeft().getY() - _bbox.getLowerLeft().getY()) /
        _cellSize);
    size_t xb = std::round(
        (qNd.bbox.getUpperRight().getX() - _bbox.getLowerLeft().getX()) /
        _cellSize);
    size_t yb = std::round(
        (qNd.bbox.getUpperRight().getY() - _bbox.getLowerLeft().getY()) /
        _cellSize);

    auto ll = getNode(xa, ya);
    if (!ll) ll = writeNd(xa, ya);

    auto lr = getNode(xb, ya);
    if (!lr) lr = writeNd(xb, ya);

    auto ul = getNode(xa, yb);
    if (!ul) ul = writeNd(xa, yb);

    auto ur = getNode(xb, yb);
    if (!ur) ur = writeNd(xb, yb);

    connectNodes(ll, ur, 1);
    connectNodes(ul, lr, 3);

    if (ll->pl().getPort(2)->getDeg() == 8) connectNodes(ll, lr, 2);
    if (ul->pl().getPort(2)->getDeg() == 8) connectNodes(ul, ur, 2);
    if (ll->pl().getPort(0)->getDeg() == 8) connectNodes(ll, ul, 0);
    if (lr->pl().getPort(0)->getDeg() == 8) connectNodes(lr, ur, 0);
  }

  prunePorts();
}

// _____________________________________________________________________________
double OctiQuadTree::ndMovePen(const CombNode* cbNd,
                               const GridNode* grNd) const {
  // the move penalty has to be at least the max cost of saving a single
  // grid hop - otherwise we could move the node closer and closer to the
  // other node without ever increasing the total cost.

  // additional penalty per grid move
  // TODO: make configurable
  double PEN = 0.5;

  // we may substitute a diagonal edge be a horizontal + 90 deg bend + vertical
  // edge
  double diagCost =
      _bendCosts[0] +
      fmin(_c.diagonalPen, _c.horizontalPen + _c.verticalPen + _bendCosts[2]);

  double vertCost =
      _bendCosts[0] +
      fmin(_c.verticalPen, _c.horizontalPen + _c.diagonalPen + _bendCosts[3]);

  double horiCost =
      _bendCosts[0] +
      fmin(_c.horizontalPen, _c.verticalPen + _c.diagonalPen + _bendCosts[3]);

  double penPerGrid = PEN + fmax(diagCost, fmax(vertCost, horiCost));

  // real distance between grid node and input node
  double d = dist(*cbNd->pl().getGeom(), *grNd->pl().getGeom());

  // distance normalized to grid length
  double gridD = d / getCellSize();

  // and multiplied per grid hop penalty
  return gridD * penPerGrid;
}
