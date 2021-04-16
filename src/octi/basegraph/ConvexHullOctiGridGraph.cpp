// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/ConvexHullOctiGridGraph.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::NodeCost;
using octi::basegraph::ConvexHullOctiGridGraph;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;

// _____________________________________________________________________________
void ConvexHullOctiGridGraph::init() {
  _ndIdx.resize(_grid.getXWidth() * _grid.getYHeight());

  // write nodes
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      if (skip(x, y)) continue;
      writeNd(x, y);
    }
  }

  // write grid edges
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      GridNode* center = getNode(x, y);
      if (!center) continue;

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
bool ConvexHullOctiGridGraph::skip(size_t x, size_t y) const {
  double xPos = _bbox.getLowerLeft().getX() + x * _cellSize;
  double yPos = _bbox.getLowerLeft().getY() + y * _cellSize;
  return !util::geo::contains({xPos, yPos}, _hull);
}

// _____________________________________________________________________________
GridNode* ConvexHullOctiGridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  auto a = _ndIdx[x * _grid.getYHeight() + y];
  if (a == 0) return 0;
  return _nds[a - 1];
}

// _____________________________________________________________________________
GridNode* ConvexHullOctiGridGraph::writeNd(size_t x, size_t y) {
  double xPos = _bbox.getLowerLeft().getX() + x * _cellSize;
  double yPos = _bbox.getLowerLeft().getY() + y * _cellSize;

  GridNode* n = addNd(DPoint(xPos, yPos));
  n->pl().setId(_nds.size());
  _nds.push_back(n);
  _ndIdx[x * _grid.getYHeight() + y] = _nds.size();
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

      if (x == 0 && (i == 5 || i == 6 || i == 7)) pen = INF;
      if (y == 0 && (i == 0 || i == 7 || i == 1)) pen = INF;
      if (x == _grid.getXWidth() - 1 && (i == 1 || i == 2 || i == 3)) pen = INF;
      if (y == _grid.getYHeight() - 1 && (i == 3 || i == 4 || i == 5))
        pen = INF;

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
