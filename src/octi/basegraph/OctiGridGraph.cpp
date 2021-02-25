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
GridNode* OctiGridGraph::neigh(size_t cx, size_t cy, size_t i) const {
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

  auto ge = getNEdg(a, b);
  auto gf = getNEdg(b, a);

  assert(ge);
  assert(gf);

  ge->pl().clearResEdges();
  gf->pl().clearResEdges();

  _resEdgs[ge] = 0;
  _resEdgs[gf] = 0;

  if (!a->pl().isSettled()) openTurns(a);
  if (!b->pl().isSettled()) openTurns(b);

  // unblock blocked diagonal edges crossing this edge
  size_t dir = getDir(a, b);
  if (dir % 2 != 0) {
    size_t x = a->pl().getX();
    size_t y = a->pl().getY();

    auto na = neigh(x, y, (dir + 7) % 8);
    auto nb = neigh(x, y, (dir + 1) % 8);

    if (na && nb) {
      auto e = getNEdg(na, nb);
      auto f = getNEdg(nb, na);

      e->pl().unblock();
      f->pl().unblock();
    }
  }
}

// _____________________________________________________________________________
size_t OctiGridGraph::ang(size_t i, size_t j) const {
  int d = (int)(i) - (int)(j);
  int deg = abs((((d + 4) % 8) + 8) % 8 - 4) % 4;

  int ang = (8 + (i - j)) % 8;
  if (ang > 4) ang = 8 - ang;
  ang = ang % 4;

  assert(deg == ang);

  return ang;
}

// _____________________________________________________________________________
double OctiGridGraph::getBendPen(size_t i, size_t j) const {
  return _bendCosts[ang(i, j)];
}

// _____________________________________________________________________________
void OctiGridGraph::settleEdg(GridNode* a, GridNode* b, CombEdge* e,
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
  size_t dir = getDir(a, b);
  if (dir % 2 != 0) {
    size_t x = a->pl().getX();
    size_t y = a->pl().getY();

    auto na = neigh(x, y, (dir + 7) % 8);
    auto nb = neigh(x, y, (dir + 1) % 8);

    if (na && nb) {
      auto e = getNEdg(na, nb);
      auto f = getNEdg(nb, na);

      e->pl().block();
      f->pl().block();
    }
  }
}

// _____________________________________________________________________________
size_t OctiGridGraph::getDir(const GridNode* a, const GridNode* b) const {
  if (!a || !b || a == b) return 0;

  int i = (int)a->pl().getX() - (int)b->pl().getX();
  if (i < -1) i = -1;
  if (i > 1) i = 1;

  int j = (int)a->pl().getY() - (int)b->pl().getY();
  if (j < -1) j = -1;
  if (j > 1) j = 1;

  int aa = 1 + i;
  int bb = 1 + j;

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

  return dir;
}

// _____________________________________________________________________________
CrossEdgPairs OctiGridGraph::getCrossEdgPairs() const {
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
void OctiGridGraph::writeInitialCosts() {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto n = getNode(x, y);
      for (size_t i = 0; i < maxDeg(); i++) {
        auto port = n->pl().getPort(i);
        auto neighbor = neigh(x, y, i);

        if (!neighbor || !port) continue;

        auto oPort = neighbor->pl().getPort((i + maxDeg() / 2) % maxDeg());
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

  size_t dir = getDir(a, b);

  if (a->pl().getPort(dir) &&
      b->pl().getPort((dir + maxDeg() / 2) % maxDeg())) {
    return const_cast<GridEdge*>(
        getEdg(a->pl().getPort(dir),
               b->pl().getPort((dir + maxDeg() / 2) % maxDeg())));
  }

  return 0;
}

// _____________________________________________________________________________
GridNode* OctiGridGraph::writeNd(size_t x, size_t y) {
  double xPos = _bbox.getLowerLeft().getX() + x * _cellSize;
  double yPos = _bbox.getLowerLeft().getY() + y * _cellSize;

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

    auto e = new GridEdge(n, nn, GridEdgePL(INF, true, true, false));
    e->pl().setId(_edgeCount);
    _edgeCount++;
    n->addEdge(e);
    nn->addEdge(e);

    e = new GridEdge(nn, n, GridEdgePL(INF, true, true, false));
    e->pl().setId(_edgeCount);
    _edgeCount++;
    n->addEdge(e);
    nn->addEdge(e);
  }

  // in-node connections
  for (size_t i = 0; i < maxDeg(); i++) {
    for (size_t j = i + 1; j < maxDeg(); j++) {
      double pen = getBendPen(i, j);

      if (x == 0 && (i == 5 || i == 6 || i == 7)) pen = INF;
      if (y == 0 && (i == 0 || i == 7 || i == 1)) pen = INF;
      if (x == _grid.getXWidth() - 1 && (i == 1 || i == 2 || i == 3)) pen = INF;
      if (y == _grid.getYHeight() - 1 && (i == 3 || i == 4 || i == 5))
        pen = INF;

      auto e = new GridEdge(n->pl().getPort(i), n->pl().getPort(j),
                            GridEdgePL(pen, true, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      e->getFrom()->addEdge(e);
      e->getTo()->addEdge(e);

      e = new GridEdge(n->pl().getPort(j), n->pl().getPort(i),
                       GridEdgePL(pen, true, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      e->getFrom()->addEdge(e);
      e->getTo()->addEdge(e);
    }
  }

  return n;
}

// _____________________________________________________________________________
size_t OctiGridGraph::maxDeg() const { return 8; }

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
  // we can take at most min(dx, dy) diagonal edges. Each diagonal edge saves us
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

// _____________________________________________________________________________
double OctiGridGraph::ndMovePen(const CombNode* cbNd,
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
