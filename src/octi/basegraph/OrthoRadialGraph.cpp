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
#include "util/geo/CircularSegment.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::NodeCost;
using octi::basegraph::OrthoRadialGraph;
using util::geo::CircularSegment;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;

// _____________________________________________________________________________
void OrthoRadialGraph::init() {
  // write nodes
  // TODO: we are only going from 1 because we have no center node
  for (size_t y = 1; y < _grid.getYHeight() / 2; y++) {
    for (size_t x = 0; x < _numBeams; x++) {
      writeNd(x, y);
    }
  }

  // write grid edges
  for (size_t x = 0; x < _numBeams; x++) {
    for (size_t y = 0; y < _grid.getYHeight() / 2; y++) {
      GridNode* center = getNode(x, y);
      if (!center) continue;

      for (size_t p = 0; p < maxDeg(); p++) {
        GridNode* from = center->pl().getPort(p);
        GridNode* toN = neigh(x, y, p);
        if (from != 0 && toN != 0) {
          GridNode* to = toN->pl().getPort((p + maxDeg() / 2) % maxDeg());
          if (!to) continue;
          auto e = new GridEdge(from, to, GridEdgePL(9, false, false));
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
GridEdge* OrthoRadialGraph::getNEdg(const GridNode* a,
                                    const GridNode* b) const {
  if (!a || !b) return 0;

  int dy = (int)a->pl().getY() - (int)b->pl().getY();
  int dx = (int)a->pl().getX() - (int)b->pl().getX();

  size_t dir = 0;

  if (dy == 1 && dx == 0)
    dir = 2;
  else if (dy == -1 && dx == 0)
    dir = 0;
  else if (dy == 0 && (dx == -((int)_numBeams - 1) || dx == 1))
    dir = 3;
  else if (dy == 0 && (dx == ((int)_numBeams - 1) || dx == -1))
    dir = 1;
  else
    return 0;

  if (a->pl().getPort(dir) &&
      b->pl().getPort((dir + maxDeg() / 2) % maxDeg())) {
    return const_cast<GridEdge*>(
        getEdg(a->pl().getPort(dir),
               b->pl().getPort((dir + maxDeg() / 2) % maxDeg())));
  }

  return 0;
}

// _____________________________________________________________________________
void OrthoRadialGraph::writeInitialCosts() {
  double angStep = 2.0 * M_PI / _numBeams;

  double c_0 = _c.p_45 - _c.p_135;

  for (size_t x = 0; x < _numBeams; x++) {
    for (size_t y = 0; y < _grid.getYHeight() / 2; y++) {
      auto n = getNode(x, y);
      if (!n) continue;
      for (size_t i = 0; i < maxDeg(); i++) {
        auto port = n->pl().getPort(i);
        auto neighbor = neigh(x, y, i);

        if (!neighbor || !port) continue;

        auto oPort = neighbor->pl().getPort((i + maxDeg() / 2) % maxDeg());
        auto e = getEdg(port, oPort);

        // this is percentage the hop length is longer than the smallest cell
        // size
        double sX = angStep / angStep * y;

        // here, the costs are normalized to always
        // represent the map lengths exactly
        if (i % 2 == 0) {
          // vertical hops always have the same length
          e->pl().setCost((_c.verticalPen + c_0) - c_0);
        } else {
          // horizontal hops get bigger with higher y (= higher radius)
          e->pl().setCost((_c.horizontalPen + c_0) * sX - c_0);
        }
      }
    }
  }
}

// _____________________________________________________________________________
double OrthoRadialGraph::ndMovePen(const CombNode* cbNd,
                                   const GridNode* grNd) const {
  // the move penalty has to be at least the max cost of saving a single
  // grid hop - otherwise we could move the node closer and closer to the
  // other node without ever increasing the total cost.

  // additional penalty per grid move
  // TODO: make configurable
  double PEN = 0.5;

  double angStep = 2.0 * M_PI / _numBeams;

  double rat = 1 / angStep;

  double c_0 = _c.p_45 - _c.p_135;
  double penPerGrid = PEN + c_0 + rat * fmax(_c.verticalPen, _c.horizontalPen);

  // real distance between grid node and input node
  double d = dist(*cbNd->pl().getGeom(), *grNd->pl().getGeom());

  // distance normalized to grid length
  double gridD = d / (getCellSize() * angStep);

  // and multiplied per grid hop penalty
  return gridD * penPerGrid;
}

// _____________________________________________________________________________
GridNode* OrthoRadialGraph::writeNd(size_t x, size_t y) {
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
  auto pos = DPoint(xPos, yPos);

  double c_0 = _c.p_45 - _c.p_135;
  double c_90 = _c.p_45 - _c.p_135 + _c.p_90;

  GridNode* n = addNd(pos);
  n->pl().setId(_nds.size());
  _nds.push_back(n);
  n->pl().setSink();
  _grid.add(pos, n);
  n->pl().setXY(x, y);
  n->pl().setParent(n);

  for (int i = 0; i < 4; i++) {
    int xi = 0;
    int yi = 0;

    if (i == 0) yi = 1;
    if (i == 1) xi = 1;
    if (i == 2) yi = -1;
    if (i == 3) xi = -1;

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

      if (y == 1 && i == 2) pen = INF;
      if (y == 1 && j == 2) pen = INF;
      if (y == _grid.getYHeight() / 2 && i == 0) pen = INF;

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
GridNode* OrthoRadialGraph::getNode(size_t x, size_t y) const {
  if (y == 0) return 0;
  if (((y - 1) * _numBeams + x) * 5 >= _nds.size()) return 0;
  return _nds[((y - 1) * _numBeams + x) * 5];
}

// _____________________________________________________________________________
PolyLine<double> OrthoRadialGraph::geomFromPath(
    const std::vector<std::pair<size_t, size_t>>& res) const {
  double origX =
      _bbox.getLowerLeft().getX() +
      (_bbox.getUpperRight().getX() - _bbox.getLowerLeft().getX()) * 0.5;

  double origY =
      _bbox.getLowerLeft().getY() +
      (_bbox.getUpperRight().getY() - _bbox.getLowerLeft().getY()) * 0.5;

  PolyLine<double> pl;
  for (auto revIt = res.rbegin(); revIt != res.rend(); revIt++) {
    auto f = getEdg(getGrNdById(revIt->first), getGrNdById(revIt->second));
    auto frPar = f->getFrom()->pl().getParent();
    auto toPar = f->getTo()->pl().getParent();

    if (!f->pl().isSecondary()) {
      pl << *frPar->pl().getGeom();

      if (frPar->pl().getY() == toPar->pl().getY()) {
        double angStep = 2.0 * M_PI / _numBeams;
        if (((frPar->pl().getX() - toPar->pl().getX()) == 1) ||
            ((toPar->pl().getX() - frPar->pl().getX()) == _numBeams - 1)) {
          angStep = -angStep;
        }

        CircularSegment<double> seg(*frPar->pl().getGeom(), -angStep,
                                    {origX, origY});
        for (auto p : seg.render(10).getLine()) pl << p;
      }

      pl << *toPar->pl().getGeom();

      pl.simplify(1);
    }
  }

  return pl;
}
