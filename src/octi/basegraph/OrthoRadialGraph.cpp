// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/OrthoRadialGraph.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::NodeCost;
using octi::basegraph::OrthoRadialGraph;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;

// _____________________________________________________________________________
void OrthoRadialGraph::init() {
  // write nodes
  for (size_t y = 0; y < _grid.getYHeight() / 2; y++) {
    for (size_t x = 0; x < _numBeams; x++) {
      writeNd(x, y);
    }
  }

  // write grid edges
  for (size_t x = 0; x < _numBeams; x++) {
    for (size_t y = 0; y < _grid.getYHeight() / 2; y++) {
      GridNode* center = getNode(x, y);

      for (size_t p = 0; p < maxDeg(); p++) {
        GridNode* from = center->pl().getPort(p);
        GridNode* toN = neigh(x, y, p);
        if (from != 0 && toN != 0) {
          GridNode* to = toN->pl().getPort((p + maxDeg() / 2) %
                                           maxDeg());
          if (!to) continue;
          auto e = new GridEdge(from, to, GridEdgePL(9, false));
          e->pl().setId(_edgeCount);
          _edgeCount++;
          from->addEdge(e);
          to->addEdge(e);
        }
      }
    }
  }

  writeInitialCosts();
}

// _____________________________________________________________________________
GridNode* OrthoRadialGraph::neigh(size_t cx, size_t cy, size_t i) const {
  if (i == 0) return getNode(cx, cy + 1);
  if (i == 1) return getNode((cx + 1) % _numBeams, cy);
  if (i == 2) return getNode(cx, cy - 1);
  if (i == 3) return getNode((cx + (_numBeams - 1)) % _numBeams, cy);

  return getNode(cx, cy);
}

// _____________________________________________________________________________
const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
OrthoRadialGraph::getHeur(const std::set<GridNode*>& to) const {
  return new OrthoRadialGraphHeur(this, to);
}

// _____________________________________________________________________________
GridEdge* OrthoRadialGraph::getNEdg(const GridNode* a, const GridNode* b) const {
  if (!a || !b) return 0;

  int dy = (int)a->pl().getY() - (int)b->pl().getY();
  int dx = (int)a->pl().getX() - (int)b->pl().getX();

  size_t dir = 0;

  if (dy == 1 && dx == 0) dir = 2;
  else if (dy == -1 && dx == 0) dir = 0;
  else if (dy == 0 && (dx == -((int)_numBeams - 1) || dx == 1)) dir = 3;
  else if (dy == 0 && (dx == ((int)_numBeams - 1) || dx == -1)) dir = 1;
  else return 0;

  // std::cerr << "(" << a->pl().getX() << "," <<  a->pl().getY()  << ")" << ", " << "(" << b->pl().getX() << "," <<  b->pl().getY()  << ")" <<  " dir: " << dir << std::endl;

  if (a->pl().getPort(dir) &&
      b->pl().getPort((dir + maxDeg() / 2) % maxDeg())) {
    return const_cast<GridEdge*>(getEdg(
        a->pl().getPort(dir),
        b->pl().getPort((dir + maxDeg() / 2) % maxDeg())));
  }

  return 0;
}

// _____________________________________________________________________________
void OrthoRadialGraph::writeInitialCosts() {
  for (size_t x = 0; x < _numBeams; x++) {
    for (size_t y = 0; y < _grid.getYHeight() / 2; y++) {
      auto n = getNode(x, y);
      for (size_t i = 0; i < maxDeg(); i++) {
        auto port = n->pl().getPort(i);
        auto neighbor = neigh(x, y, i);

        if (!neighbor || !port) continue;

        auto oPort = neighbor->pl().getPort((i + maxDeg() / 2) %
                                         maxDeg());
        auto e = getEdg(port, oPort);

        if (i % 2 == 0) {
          e->pl().setCost(_c.verticalPen);
        } else {
          e->pl().setCost(_c.horizontalPen);
        }
      }
    }
  }
}

// _____________________________________________________________________________
GridNode* OrthoRadialGraph::writeNd(size_t x, size_t y) {
  y += 1;
  double origX =
      _bbox.getLowerLeft().getX() +
      (_bbox.getUpperRight().getX() - _bbox.getLowerLeft().getX()) * 0.5;

  double origY =
      _bbox.getLowerLeft().getY() +
      (_bbox.getUpperRight().getY() - _bbox.getLowerLeft().getY()) * 0.5;

  double angStep = 2.0 * M_PI / _numBeams;

  double curAngle = -angStep * x + 0.5 * M_PI;

  double xPos = origX + y * cos(curAngle) * _cellSize;
  double yPos = origY + y * sin(curAngle) * _cellSize;

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

  for (int i = 0; i < 4; i++) {
    int xi = 0;
    int yi = 0;

    if (i == 0) {
      yi = 1;
    }
    if (i == 1) {
      xi = 1;
    }
    if (i == 2) {
      yi = -1;
    }
    if (i == 3) {
      xi = -1;
    }

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
  for (size_t i = 0; i < maxDeg(); i++) {
    for (size_t j = i + 1; j < maxDeg(); j++) {
      int d = (int)(i) - (int)(j);
      size_t deg = abs((((d + 2) % 4) + 4) % 4 - 2);
      double pen = c_0;

      if (deg == 0) pen = c_0;
      if (deg == 1) pen = c_90;

      if (x == 0 && i == 3) pen = INF;
      if (y == 0 && i == 0) pen = INF;
      if (x == _grid.getXWidth() - 1 && i == 1) pen = INF;
      if (y == _grid.getYHeight() - 1 && i == 2) pen = INF;

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
GridNode* OrthoRadialGraph::getNode(size_t x, size_t y) const {
  if ((y * _numBeams + x) * 5 >= _nds.size()) return 0;
  return _nds[(y * _numBeams + x) * 5];
}
