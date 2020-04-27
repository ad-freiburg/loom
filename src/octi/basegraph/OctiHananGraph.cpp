// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/OctiHananGraph.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::NodeCost;
using octi::basegraph::OctiHananGraph;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;

// _____________________________________________________________________________
GridNode* OctiHananGraph::neigh(size_t cx, size_t cy, size_t i) const {
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
void OctiHananGraph::unSettleEdg(GridNode* a, GridNode* b) {
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
size_t OctiHananGraph::ang(size_t i, size_t j) const {
  int d = (int)(i) - (int)(j);
  int deg = abs((((d + 4) % 8) + 8) % 8 - 4) % 4;

  int ang = (8 + (i - j)) % 8;
  if (ang > 4) ang = 8 - ang;
  ang = ang % 4;

  assert(deg == ang);

  return ang;
}

// _____________________________________________________________________________
double OctiHananGraph::getBendPen(size_t i, size_t j) const {
  return _bendCosts[ang(i, j)];
}

// _____________________________________________________________________________
void OctiHananGraph::settleEdg(GridNode* a, GridNode* b, CombEdge* e,
                              size_t rndrOrd) {
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

  ge->pl().setRndrOrder(rndrOrd);

  // this closes both nodes
  // a close means that all major edges reaching this node are closed
  closeTurns(a);
  closeTurns(b);

  // block diagonal edges crossing this edge
  if (dir != 0) {
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
CrossEdgPairs OctiHananGraph::getCrossEdgPairs() const {
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
void OctiHananGraph::writeInitialCosts() {
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
GridEdge* OctiHananGraph::getNEdg(const GridNode* a, const GridNode* b) const {
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
      b->pl().getPort((dir + maxDeg() / 2) % maxDeg())) {
    return const_cast<GridEdge*>(
        getEdg(a->pl().getPort(dir),
               b->pl().getPort((dir + maxDeg() / 2) % maxDeg())));
  }

  return 0;
}

// _____________________________________________________________________________
void OctiHananGraph::init() {

  std::vector<GridNode*> xSorted;
  std::vector<GridNode*> ySorted;

  std::map<const CombNode*, std::pair<size_t, size_t>> coords;

  // write nodes
  for (auto cNd : _cg.getNds()) {
    int x = _grid.getCellXFromX(cNd->pl().getGeom()->getX());
    int y = _grid.getCellYFromY(cNd->pl().getGeom()->getY());

    if (!getNode(x, y)) {
      auto nd = writeNd(x, y);
      coords[cNd] = {x, y};

      xSorted.push_back(nd);
      ySorted.push_back(nd);
    }
  }

  struct {
    bool operator() (const GridNode* a, const GridNode* b) {
      return a->pl().getX() < b->pl().getX();
    }
  } sortByX;

  struct {
    bool operator() (const GridNode* a, const GridNode* b) {
      return a->pl().getY() < b->pl().getY();
    }
  } sortByY;

  if (xSorted.size() == 0) return;

  std::vector<std::vector<GridNode*>> yAct(_grid.getYHeight());
  std::vector<std::vector<GridNode*>> xAct(_grid.getXWidth());

  std::vector<std::vector<GridNode*>> xyAct(_grid.getXWidth() + _grid.getYHeight());
  std::vector<std::vector<GridNode*>> yxAct(_grid.getXWidth() + _grid.getYHeight());

  for (auto nd : xSorted) {
    yAct[nd->pl().getY()].push_back(nd);
  }

  for (auto nd : xSorted) {
    xAct[nd->pl().getX()].push_back(nd);
  }

  for (auto nd : xSorted) {
    xyAct[nd->pl().getX() + (_grid.getYHeight() - 1 - nd->pl().getY())].push_back(nd);
  }

  for (auto nd : xSorted) {
    yxAct[nd->pl().getY() + nd->pl().getX()].push_back(nd);
  }

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    if (!xAct[x].size()) continue;

    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      if (!yAct[y].size()) continue;
      if (getNode(x, y)) continue;
      auto newNd = writeNd(x, y);
      yAct[y].push_back(newNd);
      xAct[x].push_back(newNd);
    }
  }

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      size_t xi = x + (_grid.getYHeight() - 1 - y);
      size_t yi = y + x;
      if ((xyAct[xi].size() && (yxAct[yi].size() || yAct[y].size() || xAct[x].size())) ||
       (yxAct[yi].size() && (xyAct[xi].size() || yAct[y].size() || xAct[x].size()))) {

        bool have = false;

        auto newNd = getNode(x, y);

        if (newNd) {
          have = true;
        }

        if (!have) newNd = writeNd(x, y);

        if (xyAct[xi].size()) {
          xyAct[xi].push_back(newNd);
        }

        if (yxAct[yi].size()) {
          yxAct[yi].push_back(newNd);
        }

        if (have) continue;

        if (yAct[y].size()) {
          yAct[y].push_back(newNd);
        }

        if (xAct[x].size()) {
          xAct[x].push_back(newNd);
        }

      }
    }
  }

  for (size_t x = 0; x < _grid.getYHeight(); x++) {
    std::sort(xAct[x].begin(), xAct[x].end(), sortByY);
  }

  for (size_t y = 0; y < _grid.getYHeight(); y++) {
    std::sort(yAct[y].begin(), yAct[y].end(), sortByX);
  }

  for (size_t i = 0; i < _grid.getYHeight() + _grid.getXWidth(); i++) {
    std::sort(xyAct[i].begin(), xyAct[i].end(), sortByY);
    std::sort(yxAct[i].begin(), yxAct[i].end(), sortByX);
  }


  for (size_t y = 0; y < _grid.getYHeight();y++) {
    for (size_t i = 1; i < yAct[y].size(); i++) {
      auto e = new GridEdge(yAct[y][i-1], yAct[y][i], GridEdgePL(9, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      yAct[y][i-1]->addEdge(e);
      yAct[y][i]->addEdge(e);

      auto f = new GridEdge(yAct[y][i], yAct[y][i-1], GridEdgePL(9, false));
      f->pl().setId(_edgeCount);
      _edgeCount++;
      yAct[y][i-1]->addEdge(f);
      yAct[y][i]->addEdge(f);
    }
  }

  for (size_t x = 0; x < _grid.getXWidth();x++) {
    for (size_t i = 1; i < xAct[x].size(); i++) {
      auto e = new GridEdge(xAct[x][i-1], xAct[x][i], GridEdgePL(9, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      xAct[x][i-1]->addEdge(e);
      xAct[x][i]->addEdge(e);

      auto f = new GridEdge(xAct[x][i], xAct[x][i-1], GridEdgePL(9, false));
      f->pl().setId(_edgeCount);
      _edgeCount++;
      xAct[x][i-1]->addEdge(f);
      xAct[x][i]->addEdge(f);
    }
  }

  for (size_t xi = 0; xi < _grid.getXWidth() + _grid.getYHeight() ;xi++) {
    for (size_t i = 1; i < xyAct[xi].size(); i++) {
      auto e = new GridEdge(xyAct[xi][i-1], xyAct[xi][i], GridEdgePL(9, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      xyAct[xi][i-1]->addEdge(e);
      xyAct[xi][i]->addEdge(e);

      auto f = new GridEdge(xyAct[xi][i], xyAct[xi][i-1], GridEdgePL(9, false));
      f->pl().setId(_edgeCount);
      _edgeCount++;
      xyAct[xi][i-1]->addEdge(f);
      xyAct[xi][i]->addEdge(f);
    }
  }

  for (size_t yi = 0; yi < _grid.getXWidth() + _grid.getYHeight() ;yi++) {
    for (size_t i = 1; i < yxAct[yi].size(); i++) {
      auto e = new GridEdge(yxAct[yi][i-1], yxAct[yi][i], GridEdgePL(9, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;
      yxAct[yi][i-1]->addEdge(e);
      yxAct[yi][i]->addEdge(e);

      auto f = new GridEdge(yxAct[yi][i], yxAct[yi][i-1], GridEdgePL(9, false));
      f->pl().setId(_edgeCount);
      _edgeCount++;
      yxAct[yi][i-1]->addEdge(f);
      yxAct[yi][i]->addEdge(f);
    }
  }



  // write grid edges
  // for (size_t x = 0; x < _grid.getXWidth(); x++) {
    // for (size_t y = 0; y < _grid.getYHeight(); y++) {
      // GridNode* center = getNode(x, y);

      // for (size_t p = 0; p < maxDeg(); p++) {
        // GridNode* frN = center->pl().getPort(p);
        // GridNode* toN = neigh(x, y, p);
        // if (frN != 0 && toN != 0) {
          // GridNode* to = toN->pl().getPort((p + maxDeg() / 2) % maxDeg());
          // auto e = new GridEdge(frN, to, GridEdgePL(9, false));
          // e->pl().setId(_edgeCount);
          // _edgeCount++;
          // frN->addEdge(e);
          // to->addEdge(e);
        // }
      // }
    // }
  // }

  // writeInitialCosts();
}

// _____________________________________________________________________________
GridNode* OctiHananGraph::writeNd(size_t x, size_t y) {
  double xPos = _bbox.getLowerLeft().getX() + x * _cellSize;
  double yPos = _bbox.getLowerLeft().getY() + y * _cellSize;

  GridNode* n = addNd(DPoint(xPos, yPos));
  n->pl().setId(_nds.size());
  _nds[{x, y}] = n;
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
  for (size_t i = 0; i < maxDeg(); i++) {
    for (size_t j = i + 1; j < maxDeg(); j++) {
      double pen = getBendPen(i, j);

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
size_t OctiHananGraph::maxDeg() const { return 8; }

// _____________________________________________________________________________
GridNode* OctiHananGraph::getNode(size_t x, size_t y) const {
  auto f = _nds.find({x, y});
  if (f == _nds.end()) return 0;
  return f->second;
}

// _____________________________________________________________________________
double OctiHananGraph::heurCost(int64_t xa, int64_t ya, int64_t xb,
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

// _____________________________________________________________________________
size_t OctiHananGraph::getGrNdDeg(const CombNode* nd, size_t x, size_t y) const {
  if ((x == 0 || x == _grid.getXWidth() - 1) &&
      (y == 0 || y == _grid.getYHeight() - 1))
    return 3;

  if ((x == 0 || x == _grid.getXWidth() - 1) ||
      (y == 0 || y == _grid.getYHeight() - 1))
    return 5;

  return 8;
}

// _____________________________________________________________________________
double OctiHananGraph::ndMovePen(const CombNode* cbNd,
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
