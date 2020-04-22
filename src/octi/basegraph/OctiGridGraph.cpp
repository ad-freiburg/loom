// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/OctiGridGraph.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::NodeCost;
using octi::basegraph::OctiGridGraph;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;

// _____________________________________________________________________________
GridNode* OctiGridGraph::getNeighbor(size_t cx, size_t cy, size_t i) const {
  if (i > 7) return getNode(cx, cy);
  int8_t x = 1;
  if (i % 4 == 0) x = 0;
  if (i > 4) x = -1;

  int8_t y = 1;
  if (i == 2 || i == 6) y = 0;
  if (i == 3 || i == 4 || i == 5) y = -1;

  return getNode(cx + x, cy + y);
}

// _____________________________________________________________________________
void OctiGridGraph::unSettleEdg(GridNode* a, GridNode* b) {
  if (a == b) return;
  int aa = 1 + (int)a->pl().getX() - (int)b->pl().getX();
  int bb = 1 + (int)a->pl().getY() - (int)b->pl().getY();

  size_t dir = 0;
  size_t d = aa * 3 + bb;

  if (d == 0) dir = 1;
  if (d == 2) dir = 3;
  if (d == 8) dir = 5;
  if (d == 6) dir = 7;

  size_t x = a->pl().getX();
  size_t y = a->pl().getY();

  auto ge = getNEdg(a, b);
  auto gf = getNEdg(b, a);

  assert(ge);
  assert(gf);

  ge->pl().clearResEdges();
  gf->pl().clearResEdges();

  _resEdgs[ge] = 0;
  _resEdgs[gf] = 0;

  if (!a->pl().isSettled()) {
    openTurns(a);
  }
  if (!b->pl().isSettled()) {
    openTurns(b);
  }

  // unblock blocked diagonal edges crossing this edge
  if (dir != 0) {
    auto na = getNeighbor(x, y, (dir + 7) % 8);
    auto nb = getNeighbor(x, y, (dir + 1) % 8);

    if (na && nb) {
      auto e = getNEdg(na, nb);
      auto f = getNEdg(nb, na);

      e->pl().unblock();
      f->pl().unblock();
    }
  }
}

// _____________________________________________________________________________
double OctiGridGraph::getBendPen(size_t origI, size_t targetI) const {
  // determine angle between port i and j
  int ang = (8 + (origI - targetI)) % 8;
  if (ang > 4) ang = 8 - ang;
  ang = ang % 4;

  return _bendCosts[ang];
}

// _____________________________________________________________________________
void OctiGridGraph::settleEdg(GridNode* a, GridNode* b, CombEdge* e) {
  if (a == b) return;

  int aa = 1 + (int)a->pl().getX() - (int)b->pl().getX();
  int bb = 1 + (int)a->pl().getY() - (int)b->pl().getY();

  size_t dir = 0;
  size_t d = aa * 3 + bb;

  if (d == 0) dir = 1;
  if (d == 2) dir = 3;
  if (d == 8) dir = 5;
  if (d == 6) dir = 7;

  size_t x = a->pl().getX();
  size_t y = a->pl().getY();

  // this closes the grid edge
  auto ge = getNEdg(a, b);

  addResEdg(ge, e);

  // this closes both nodes
  // a close means that all major edges reaching this node are closed
  closeTurns(a);
  closeTurns(b);

  // block diagonal edges crossing this edge
  if (dir != 0) {
    auto na = getNeighbor(x, y, (dir + 7) % 8);
    auto nb = getNeighbor(x, y, (dir + 1) % 8);

    if (na && nb) {
      auto e = getNEdg(na, nb);
      auto f = getNEdg(nb, na);

      e->pl().block();
      f->pl().block();
    }
  }
}

// _____________________________________________________________________________
CrossEdgPairs OctiGridGraph::getCrossEdgPairs() const {
  CrossEdgPairs ret;
  for (const GridNode* n : getNds()) {
    if (!n->pl().isSink()) continue;

    auto eOr = getNEdg(n, getNeighbor(n, 3));
    auto fOr = getNEdg(getNeighbor(n, 3), n);

    if (!eOr || !fOr) continue;

    auto na = getNeighbor(n, (3 + 7) % 8);
    auto nb = getNeighbor(n, (3 + 1) % 8);

    if (!na || !nb) continue;

    auto e = getNEdg(na, nb);
    auto f = getNEdg(nb, na);
    ret.push_back({{eOr, fOr}, {e, f}});
  }

  return ret;
}

// _____________________________________________________________________________
void OctiGridGraph::writeInitialCosts() {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto n = getNode(x, y);
      for (size_t i = 0; i < getNumNeighbors(); i++) {
        auto port = n->pl().getPort(i);
        auto neigh = getNeighbor(x, y, i);

        if (!neigh || !port) continue;

        auto oPort = neigh->pl().getPort((i + getNumNeighbors() / 2) %
                                         getNumNeighbors());
        auto e = getEdg(port, oPort);

        if (i % 4 == 0) {
          e->pl().setCost(_c.verticalPen);
        } else if ((i + 2) % 4 == 0) {
          e->pl().setCost(_c.horizontalPen);
        } else if (i % 2) {
          e->pl().setCost(_c.diagonalPen);
        }
      }
    }
  }
}

// _____________________________________________________________________________
GridEdge* OctiGridGraph::getNEdg(const GridNode* a, const GridNode* b) const {
  if (!a || !b) return 0;

  int aa = 1 + (int)a->pl().getX() - (int)b->pl().getX();
  int bb = 1 + (int)a->pl().getY() - (int)b->pl().getY();

  size_t dir = 0;
  size_t d = aa * 3 + bb;

  if (d == 0) dir = 1;
  if (d == 1) dir = 2;
  if (d == 2) dir = 3;
  if (d == 3) dir = 0;
  if (d == 5) dir = 4;
  if (d == 8) dir = 5;
  if (d == 7) dir = 6;
  if (d == 6) dir = 7;

  if (a->pl().getPort(dir) &&
      b->pl().getPort((dir + getNumNeighbors() / 2) % getNumNeighbors())) {
    return const_cast<GridEdge*>(getEdg(
        a->pl().getPort(dir),
        b->pl().getPort((dir + getNumNeighbors() / 2) % getNumNeighbors())));
  }

  return 0;
}

// _____________________________________________________________________________
GridNode* OctiGridGraph::writeNd(size_t x, size_t y) {
  double xPos = _bbox.getLowerLeft().getX() + x * _cellSize;
  double yPos = _bbox.getLowerLeft().getY() + y * _cellSize;

  double c_0 = _c.p_45 - _c.p_135;
  double c_135 = _c.p_45;
  double c_90 = _c.p_45 - _c.p_135 + _c.p_90;
  double c_45 = c_0 + c_135;

  GridNode* n = addNd(DPoint(xPos, yPos));
  n->pl().setId(_nds.size());
  _nds.push_back(n);
  n->pl().setSink();
  _grid.add(x, y, n);
  n->pl().setXY(x, y);
  n->pl().setParent(n);

  for (int i = 0; i < 8; i++) {
    int xi = (4 - (i % 8)) % 4;
    xi /= abs(abs(xi) - 1) + 1;
    int yi = (4 - ((i + 2) % 8)) % 4;
    yi /= abs(abs(yi) - 1) + 1;

    GridNode* nn = addNd(DPoint(xPos + xi * _spacer, yPos + yi * _spacer));
    nn->pl().setId(_nds.size());
    _nds.push_back(nn);
    nn->pl().setParent(n);
    n->pl().setPort(i, nn);

    auto e = new GridEdge(n, nn, GridEdgePL(INF, true, false));
    e->pl().setId(_edgeCount);
    _edgeCount++;
    n->addEdge(e);
    nn->addEdge(e);

    e = new GridEdge(nn, n, GridEdgePL(INF, true, false));
    e->pl().setId(_edgeCount);
    _edgeCount++;
    n->addEdge(e);
    nn->addEdge(e);
  }

  // in-node connections
  for (size_t i = 0; i < getNumNeighbors(); i++) {
    for (size_t j = i + 1; j < getNumNeighbors(); j++) {
      int d = (int)(i) - (int)(j);
      size_t deg = abs((((d + 4) % 8) + 8) % 8 - 4) % 4;
      double pen = c_0;

      if (deg == 0) pen = c_0;
      if (deg == 1) pen = c_45;
      if (deg == 2) pen = c_90;
      if (deg == 3) pen = c_135;

      std::cerr << deg << " " << pen << " VS " << _bendCosts[deg] << std::endl;
      assert(pen == _bendCosts[deg]);

      // pen = _bendCosts[deg];

      if (x == 0 && (i == 5 || i == 6 || i == 7)) pen = INF;
      if (y == 0 && (i == 0 || i == 7 || i == 1)) pen = INF;
      if (x == _grid.getXWidth() - 1 && (i == 1 || i == 2 || i == 3)) pen = INF;
      if (y == _grid.getYHeight() - 1 && (i == 3 || i == 4 || i == 5))
        pen = INF;

      auto e = new GridEdge(n->pl().getPort(i), n->pl().getPort(j),
                            GridEdgePL(pen, true));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      e->getFrom()->addEdge(e);
      e->getTo()->addEdge(e);

      e = new GridEdge(n->pl().getPort(j), n->pl().getPort(i),
                       GridEdgePL(pen, true));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      e->getFrom()->addEdge(e);
      e->getTo()->addEdge(e);
    }
  }

  return n;
}

// _____________________________________________________________________________
size_t OctiGridGraph::getNumNeighbors() const { return 8; }

// _____________________________________________________________________________
GridNode* OctiGridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  return _nds[_grid.getYHeight() * 9 * x + y * 9];
}

// _____________________________________________________________________________
double OctiGridGraph::heurCost(int64_t xa, int64_t ya, int64_t xb,
                               int64_t yb) const {
  int dx = labs(xb - xa);
  int dy = labs(yb - ya);

  // cost without using diagonals
  // we can take at most min(dx, dy) diagonal edges. Edge diagonal edge saves us
  // one horizontal and one vertical edge, but costs a diagonal edge
  double edgeCost =
      _heurXCost * dx + _heurYCost * dy + _heurDiagSave * std::min(dx, dy);

  // we have to do at least one turn!
  if (dx != dy && dx != 0 && dy != 0) edgeCost += _c.p_135;

  // // Worse alternative: use a chebyshev distance heuristic
  // double minHops = std::max(dx, dy);
  // double heurECost =
  // (std::min(_c.verticalPen, std::min(_c.horizontalPen, _c.diagonalPen)));

  // double cc = minHops * (heurECost + _heurHopCost) - _heurHopCost;

  // return cc;

  // we always count one heurHopCost too much, subtract it at the end!
  return edgeCost - _heurHopCost;
}
