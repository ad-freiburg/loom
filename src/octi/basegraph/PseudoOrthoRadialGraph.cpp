// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include "octi/basegraph/NodeCost.h"
#include "octi/basegraph/PseudoOrthoRadialGraph.h"
#include "util/Misc.h"
#include "util/geo/CircularSegment.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::NodeCost;
using octi::basegraph::PseudoOrthoRadialGraph;
using util::geo::CircularSegment;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;

// _____________________________________________________________________________
int PseudoOrthoRadialGraph::multi(size_t y) const {
  int lg = log2(y);
  return 1 << lg;
}

// _____________________________________________________________________________
void PseudoOrthoRadialGraph::writeObstacleCost(
    const util::geo::Polygon<double>& obst) {
  for (size_t y = 1; y < _grid.getYHeight() / 2; y++) {
    for (size_t x = 0; x < _numBeams * multi(y); x++) {
      auto grNdA = getNode(x, y);

      for (size_t i = 0; i < maxDeg(); i++) {
        auto grNeigh = neigh(x, y, i);
        if (!grNeigh) continue;
        auto ge = getNEdg(grNdA, grNeigh);

        if (intersects(
                util::geo::LineSegment<double>(*ge->getFrom()->pl().getGeom(),
                                               *ge->getTo()->pl().getGeom()),
                obst) ||
            contains(
                util::geo::LineSegment<double>(*ge->getFrom()->pl().getGeom(),
                                               *ge->getTo()->pl().getGeom()),
                obst)) {
          ge->pl().setCost(std::numeric_limits<double>::infinity());
        }
      }
    }
  }
}

// _____________________________________________________________________________
void PseudoOrthoRadialGraph::writeGeoCoursePens(const CombEdge* ce,
                                                GeoPensMap* target,
                                                double pen) {
  std::set<GridNode*> neighs;

  DBox box;

  std::vector<util::geo::DLine> geoms;

  for (auto orE : ce->pl().getChilds()) {
    box = util::geo::extendBox(*orE->pl().getGeom(), box);

    // operate on simplified geometries
    geoms.push_back(util::geo::simplify(*orE->pl().getGeom(), 5));
  }

  box = util::geo::pad(box, sqrt(SOFT_INF / pen) * getCellSize());
  _grid.get(box, &neighs);

  for (auto grNdA : neighs) {
    for (size_t i = 0; i < maxDeg(); i++) {
      auto grNeigh = neigh(grNdA->pl().getX(), grNdA->pl().getY(), i);
      if (!grNeigh) continue;
      auto ge = getNEdg(grNdA, grNeigh);

      double d = std::numeric_limits<double>::infinity();

      for (const auto& geom : geoms) {
        double dLoc = fmax(
            fmax(dist(geom, *ge->getFrom()->pl().getGeom()) / getCellSize(),
                 dist(geom, util::geo::centroid(util::geo::MultiPoint<double>{
                                *ge->getFrom()->pl().getGeom(),
                                *ge->getFrom()->pl().getGeom()})) /
                     getCellSize()),
            dist(geom, *ge->getTo()->pl().getGeom()) / getCellSize());

        if (dLoc < d) d = dLoc;
      }

      d *= pen * d;

      if (d <= SOFT_INF) (*target)[ce][ge->pl().getId()] = d;
    }
  }
}

// _____________________________________________________________________________
void PseudoOrthoRadialGraph::init() {
  // write nodes
  for (size_t y = 1; y < _grid.getYHeight() / 2; y++) {
    for (size_t x = 0; x < _numBeams * multi(y); x++) {
      writeNd(x, y);
    }
  }

  // write center node
  writeNd(0, 0);

  // write grid edges
  for (size_t y = 0; y < _grid.getYHeight() / 2; y++) {
    size_t n = _numBeams * multi(y);
    if (y == 0) n = 1;
    for (size_t x = 0; x < n; x++) {
      GridNode* grNd = getNode(x, y);
      if (!grNd) continue;

      for (size_t p = 0; p < maxDeg(); p++) {
        GridNode* from = grNd->pl().getPort(p);
        GridNode* toN = neigh(x, y, p);
        if (from != 0 && toN != 0) {
          GridNode* to = toN->pl().getPort((p + maxDeg() / 2) % maxDeg());
          if (y == 0) to = toN->pl().getPort(2);
          if (toN->pl().getY() == 0) to = toN->pl().getPort(x / 2);
          if (!to) continue;
          auto e = addEdg(from, to, GridEdgePL(9, false, false));
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
void PseudoOrthoRadialGraph::getSettledAdjEdgs(GridNode* n, CombNode* origNd,
                                               CombEdge* outgoing[8]) {
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t i = 0; i < maxDeg(); i++) {
    outgoing[i] = 0;
    auto p = n->pl().getPort(i);
    if (!p) continue;
    auto neighbor = neigh(x, y, i);
    if (!neighbor) continue;

    auto neighP = neighbor->pl().getPort((i + maxDeg() / 2) % maxDeg());

    // special handling for center node
    if (y == 0) neighP = neighbor->pl().getPort(2);
    if (neighbor->pl().getY() == 0) neighP = neighbor->pl().getPort(x / 2);

    assert(neighP);

    auto e = getEdg(p, neighP);
    auto f = getEdg(neighP, p);
    assert(e);
    assert(f);
    auto resEdgs = getResEdgs(e);
    if (!resEdgs.size()) resEdgs = getResEdgs(f);

    if (resEdgs.size()) {
      for (auto e : resEdgs) {
        // they may be incorrect resident edges because of relaxed constraints
        if (e->getFrom() == origNd || e->getTo() == origNd) {
          outgoing[i] = e;
          break;
        }
      }
    }
  }
}

// _____________________________________________________________________________
GridNode* PseudoOrthoRadialGraph::neigh(size_t cx, size_t cy, size_t i) const {
  if (cx == 0 && cy == 0) return getNode(i * 2, 1);
  if (i == 0) {
    if (multi(cy) != multi(cy + 1)) {
      return getNode(cx * 2, cy + 1);
    }
    return getNode(cx, cy + 1);
  }
  if (i == 1) return getNode((cx + 1) % (_numBeams * multi(cy)), cy);
  if (i == 2) {
    if (cy == 1) {
      if (cx % 2)
        return 0;
      else
        return getNode(0, 0);
    }
    if (multi(cy) != multi(cy - 1)) {
      if (cx % 2) return 0;
      return getNode(cx / 2, cy - 1);
    }
    return getNode(cx, cy - 1);
  }
  if (i == 3)
    return getNode((cx + (_numBeams * multi(cy) - 1)) % (_numBeams * multi(cy)),
                   cy);

  return getNode(cx, cy);
}

// _____________________________________________________________________________
const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
PseudoOrthoRadialGraph::getHeur(const std::set<GridNode*>& to) const {
  return new PseudoOrthoRadialGraphHeur(this, to);
}

// _____________________________________________________________________________
GridEdge* PseudoOrthoRadialGraph::getNEdg(const GridNode* a,
                                          const GridNode* b) const {
  if (!a || !b) return 0;

  int dy = (int)a->pl().getY() - (int)b->pl().getY();
  int dx = (int)a->pl().getX() - (int)b->pl().getX();

  if (a->pl().getX() == 0 && a->pl().getY() == 0) {
    assert(b->pl().getY() == 1);
    assert(b->pl().getX() % 2 == 0);
    assert(getEdg(a->pl().getPort(b->pl().getX() / 2), b->pl().getPort(2)));
    return const_cast<GridEdge*>(
        getEdg(a->pl().getPort(b->pl().getX() / 2), b->pl().getPort(2)));
  }

  if (b->pl().getX() == 0 && b->pl().getY() == 0) {
    assert(a->pl().getY() == 1);
    assert(a->pl().getX() % 2 == 0);
    assert(getEdg(a->pl().getPort(2), b->pl().getPort(a->pl().getX() / 2)));
    return const_cast<GridEdge*>(
        getEdg(a->pl().getPort(2), b->pl().getPort(a->pl().getX() / 2)));
  }

  size_t dir = 0;

  if (dy == 0) {
    if (dx == -((int)_numBeams * multi(a->pl().getY()) - 1) || dx == 1)
      dir = 3;
    else if (dx == ((int)_numBeams * multi(a->pl().getY()) - 1) || dx == -1)
      dir = 1;
  } else if (dy == 1) {
    dir = 2;
  } else if (dy == -1) {
    dir = 0;
  } else {
    return 0;
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
void PseudoOrthoRadialGraph::writeInitialCosts() {
  // angle step at inner ring
  double angStep = 2.0 * M_PI / _numBeams;

  double c_0 = _c.p_45 - _c.p_135;

  for (size_t y = 0; y < _grid.getYHeight() / 2; y++) {
    size_t n = _numBeams * multi(y);
    if (y == 0) n = 1;
    double angStepLoc = 2.0 * M_PI / n;
    for (size_t x = 0; x < n; x++) {
      auto n = getNode(x, y);
      if (!n) continue;
      for (size_t i = 0; i < maxDeg(); i++) {
        auto port = n->pl().getPort(i);
        if (!port) continue;
        auto neighbor = neigh(x, y, i);
        if (!neighbor) continue;

        auto oPort = neighbor->pl().getPort((i + maxDeg() / 2) % maxDeg());
        if (y == 0) oPort = neighbor->pl().getPort(2);
        if (neighbor->pl().getY() == 0) oPort = neighbor->pl().getPort(x / 2);
        auto e = getEdg(port, oPort);

        // this is the percentage the hop length is longer than the smallest
        // cell size
        double sX = (angStepLoc * y) / angStep;
        if (y == 0) sX = 1;

        // here, the costs are normalized to always
        // represent the map lengths exactly
        if (i % 2 == 0) {
          // vertical hops always have the same length
          e->pl().setCost(_c.verticalPen);
          assert(e->pl().cost() >= 0);
        } else {
          // horizontal hops get bigger with higher y (= higher radius)
          e->pl().setCost((_c.horizontalPen + c_0) * sX - c_0);
          assert(e->pl().cost() >= 0);
        }
      }
    }
  }
}

// _____________________________________________________________________________
double PseudoOrthoRadialGraph::ndMovePen(const CombNode* cbNd,
                                         const GridNode* grNd) const {
  // the move penalty has to be at least the max cost of saving a single
  // grid hop - otherwise we could move the node closer and closer to the
  // other node without ever increasing the total cost.

  // additional penalty per grid move
  // TODO: make configurable
  double PEN = _c.ndMovePen;

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
GridNode* PseudoOrthoRadialGraph::writeNd(size_t x, size_t y) {
  double origX =
      _bbox.getLowerLeft().getX() +
      (_bbox.getUpperRight().getX() - _bbox.getLowerLeft().getX()) * 0.5;

  double origY =
      _bbox.getLowerLeft().getY() +
      (_bbox.getUpperRight().getY() - _bbox.getLowerLeft().getY()) * 0.5;

  double angStep = 2.0 * M_PI / (_numBeams * multi(y));

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

  // we are using the raw position here, as grid cells do not reflect the
  // positions in the grid graph as in the octilinear case
  _grid.add(pos, n);
  n->pl().setXY(x, y); n->pl().setParent(n);

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

    auto e = addEdg(n, nn, GridEdgePL(INF, true, false));
    e->pl().setId(_edgeCount);
    _edgeCount++;

    e = addEdg(nn, n, GridEdgePL(INF, true, false));
    e->pl().setId(_edgeCount);
    _edgeCount++;
  }

  // in-node connections
  for (size_t i = 0; i < maxDeg(); i++) {
    for (size_t j = i + 1; j < maxDeg(); j++) {
      int d = (int)(i) - (int)(j);
      size_t deg = abs((((d + 2) % 4) + 4) % 4 - 2);
      double pen = c_0;

      if (deg == 0) pen = c_0;
      if (deg == 1) pen = c_90;

      if (y == 1 && x % 2 && i == 2) pen = INF;
      if (y == 1 && x % 2 && j == 2) pen = INF;
      if (y == _grid.getYHeight() / 2 && i == 0) pen = INF;

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
GridNode* PseudoOrthoRadialGraph::getNode(size_t x, size_t y) const {
  if (x == 0 && y == 0) return _nds[_nds.size() - 5];
  int a = 0;

  for (size_t i = 1; i < y; i++) a += _numBeams * multi(i);

  if (y == 0) return 0;
  if ((a + x) * 5 >= _nds.size() - 5) return 0;
  return _nds[(a + x) * 5];
}

// _____________________________________________________________________________
PolyLine<double> PseudoOrthoRadialGraph::geomFromPath(
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
        double angStep = 2.0 * M_PI / (_numBeams * multi(frPar->pl().getY()));
        if (((frPar->pl().getX() - toPar->pl().getX()) == 1) ||
            ((toPar->pl().getX() - frPar->pl().getX()) ==
             (_numBeams * multi(frPar->pl().getY())) - 1)) {
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

// _____________________________________________________________________________
double PseudoOrthoRadialGraph::heurCost(int64_t xa, int64_t ya, int64_t xb,
                                        int64_t yb) const {
  UNUSED(xa);
  UNUSED(xb);
  int dy = labs(yb - ya);

  double edgCost = (_c.verticalPen + _heurHopCost) * dy;

  // we always count one heurHopCost too much, subtract it at the end, but
  // dont make negative
  return fmax(0, edgCost - _heurHopCost);
}
