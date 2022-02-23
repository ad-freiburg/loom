// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/HexGridGraph.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::HexGridGraph;
using util::geo::DBox;
using util::geo::DPoint;

// _____________________________________________________________________________
void HexGridGraph::init() {
  // write nodes
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      writeNd(x, y);
    }
  }

  // write grid edges
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      GridNode* center = getNode(x, y);

      for (size_t p = 0; p < maxDeg(); p++) {
        GridNode* frN = center->pl().getPort(p);
        GridNode* toN = neigh(x, y, p);
        if (frN && toN) {
          GridNode* to = toN->pl().getPort((p + maxDeg() / 2) % maxDeg());
          auto e = addEdg(frN, to, GridEdgePL(9, false, false));
          e->pl().setId(_edgeCount);
          _edgeCount++;
        }
      }
    }
  }

  writeInitialCosts();
  prunePorts();
}

// _____________________________________________________________________________
GridNode* HexGridGraph::neigh(size_t cx, size_t cy, size_t i) const {
  if (i > 5) return getNode(cx, cy);
  int x = 0;
  int y = 0;

  if (!(cy % 2)) {
    if (i == 0) {
      y = 1;
    }
    if (i == 1) {
      x = 1;
    }
    if (i == 2) {
      y = -1;
    }
    if (i == 3) {
      y = -1;
      x = -1;
    }
    if (i == 4) {
      x = -1;
    }
    if (i == 5) {
      x = -1;
      y = 1;
    }
  } else {
    if (i == 0) {
      y = 1;
      x = 1;
    }
    if (i == 1) {
      x = 1;
    }
    if (i == 2) {
      y = -1;
      x = 1;
    }
    if (i == 3) {
      y = -1;
    }
    if (i == 4) {
      x = -1;
    }
    if (i == 5) {
      y = 1;
    }
  }

  return getNode(cx + x, cy + y);
}

// _____________________________________________________________________________
size_t HexGridGraph::maxDeg() const { return 6; }

// _____________________________________________________________________________
const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
HexGridGraph::getHeur(const std::set<GridNode*>& to) const {
  return new HexGridGraphHeur(this, to);
}

// _____________________________________________________________________________
GridEdge* HexGridGraph::getNEdg(const GridNode* a,
                                    const GridNode* b) const {
  if (!a || !b) return 0;

  int aa = (int)a->pl().getX() - (int)b->pl().getX();
  int bb = (int)a->pl().getY() - (int)b->pl().getY();

  size_t dir = 0;

  if (bb == 0) {
    if (aa == -1) dir = 1;
    if (aa == 1) dir = 4;
  } else if (aa == 0) {
    if (!(a->pl().getY() % 2)) {
      if (bb == -1) dir = 0;
      if (bb == 1) dir = 2;
    } else {
      if (bb == -1) dir = 5;
      if (bb == 1) dir = 3;
    }
  } else {
    if (!(a->pl().getY() % 2)) {
      if (bb == -1 && aa == 1) dir = 5;
      if (bb == 1 && aa == 1) dir = 3;
    } else {
      if (bb == -1 && aa == 1) dir = 0;
      if (bb == 1 && aa == -1) dir = 2;
    }
  }

  if (a->pl().getPort(dir) &&
      b->pl().getPort((dir + maxDeg() / 2) % maxDeg())) {
    return const_cast<GridEdge*>(
        getEdg(a->pl().getPort(dir),
               b->pl().getPort((dir + maxDeg() / 2) % maxDeg())));
  }

  return 0;
}

// _____________________________________________________________________________
void HexGridGraph::writeInitialCosts() {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto n = getNode(x, y);
      for (size_t i = 0; i < maxDeg(); i++) {
        auto port = n->pl().getPort(i);
        auto neighbor = neigh(x, y, i);

        if (!neighbor || !port) continue;

        auto oPort = neighbor->pl().getPort((i + maxDeg() / 2) % maxDeg());
        auto e = getEdg(port, oPort);

        if (i == 1 || i == 4) {
          e->pl().setCost(_c.horizontalPen);
        } else {
          e->pl().setCost(_c.diagonalPen);
        }
      }
    }
  }
}

// _____________________________________________________________________________
GridNode* HexGridGraph::writeNd(size_t x, size_t y) {
  double xPos = _bbox.getLowerLeft().getX() + x * _a;
  double yPos = _bbox.getLowerLeft().getY() + y * _h;

  if (y % 2) xPos += _a / 2;

  auto pos = DPoint(xPos, yPos);
  GridNode* n = addNd(pos);
  n->pl().setId(_nds.size());
  _nds.push_back(n);
  n->pl().setSink();

  // we are using the raw position here, as grid cells do not reflect the
  // positions in the grid graph as in the octilinear case
  _grid.add(pos, n);
  n->pl().setXY(x, y);
  n->pl().setParent(n);

  for (int i = 0; i < 6; i++) {
    double xi = 0;
    double yi = 0;

    if (i == 0) {
      xi = .5;
      yi = A;
    }
    if (i == 1) {
      xi = A;
    }
    if (i == 2) {
      xi = .5;
      yi = -A;
    }
    if (i == 3) {
      xi = -.5;
      yi = -A;
    }
    if (i == 4) {
      xi = -1;
    }
    if (i == 5) {
      xi = -.5;
      yi = A;
    }

    GridNode* nn = addNd(DPoint(xPos + xi * _spacer, yPos + yi * _spacer));
    nn->pl().setId(_nds.size());
    _nds.push_back(nn);
    nn->pl().setParent(n);
    n->pl().setPort(i, nn);

    auto e = addEdg(n, nn, GridEdgePL(INF, true, true));
    e->pl().setId(_edgeCount);
    _edgeCount++;

    e = addEdg(nn, n, GridEdgePL(INF, true, true));
    e->pl().setId(_edgeCount);
    _edgeCount++;
  }

  // in-node connections
  for (size_t i = 0; i < maxDeg(); i++) {
    for (size_t j = i + 1; j < maxDeg(); j++) {
      double pen = getBendPen(i, j);

      if (x == 0 && i == 3) pen = INF;
      if (y == 0 && i == 0) pen = INF;
      if (x == _grid.getXWidth() - 1 && i == 1) pen = INF;
      if (y == _grid.getYHeight() - 1 && i == 2) pen = INF;

      auto e = addEdg(n->pl().getPort(i), n->pl().getPort(j),
                            GridEdgePL(pen, true, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;

      e = addEdg(n->pl().getPort(j), n->pl().getPort(i),
                       GridEdgePL(pen, true, false));
      e->pl().setId(_edgeCount);
      _edgeCount++;
    }
  }

  return n;
}

// _____________________________________________________________________________
std::vector<double> HexGridGraph::getCosts() const {
  std::vector<double> ret(3);
  ret[0] = _bendCosts[0];
  ret[1] = _bendCosts[2];
  ret[2] = _bendCosts[1];
  return ret;
 }

// _____________________________________________________________________________
double HexGridGraph::getBendPen(size_t i, size_t j) const {
  return _bendCosts[ang(i, j)];
}

// _____________________________________________________________________________
size_t HexGridGraph::ang(size_t i, size_t j) const {
  // determine angle between port i and j
  int ang = (6 + (i - j)) % 6;
  if (ang > 3) ang = 6 - ang;
  ang = ang % 3;

  return ang;
}

// _____________________________________________________________________________
GridNode* HexGridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  return _nds[_grid.getYHeight() * 7 * x + y * 7];
}
