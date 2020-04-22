// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/basegraph/GridGraph.h"
#include "octi/basegraph/NodeCost.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::basegraph;
using octi::basegraph::GridGraph;
using octi::basegraph::NodeCost;
using util::geo::contains;
using util::geo::DBox;
using util::geo::dist;
using util::geo::DPoint;
using util::geo::intersects;
using util::geo::LineSegment;

// _____________________________________________________________________________
GridGraph::GridGraph(const DBox& bbox, double cellSize, double spacer,
                     const Penalties& pens)
    : _bbox(bbox),
      _c(pens),
      _grid(cellSize, cellSize, bbox, false),
      _cellSize(cellSize),
      _spacer(spacer),
      _edgeCount(0) {
  assert(_c.p_0 <= _c.p_135);
  assert(_c.p_135 <= _c.p_90);
  assert(_c.p_90 <= _c.p_45);

  _bendCosts[0] = _c.p_45 - _c.p_135;
  _bendCosts[1] = _c.p_45 - _c.p_135 + _c.p_90;

  // cut off illegal spacer values
  if (spacer > cellSize / 2) _spacer = cellSize / 2;

  _heurHopCost = _c.p_45 - _c.p_135;
}

// _____________________________________________________________________________
void GridGraph::init() {
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
        if (frN != 0 && toN != 0) {
          GridNode* to = toN->pl().getPort((p + maxDeg() / 2) % maxDeg());
          auto e = new GridEdge(frN, to, GridEdgePL(9, false));
          e->pl().setId(_edgeCount);
          _edgeCount++;
          frN->addEdge(e);
          to->addEdge(e);
        }
      }
    }
  }

  writeInitialCosts();
}

// _____________________________________________________________________________
size_t GridGraph::maxDeg() const { return 4; }

// _____________________________________________________________________________
GridNode* GridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  return _nds[_grid.getYHeight() * 5 * x + y * 5];
}

// _____________________________________________________________________________
GridNode* GridGraph::neigh(size_t cx, size_t cy, size_t i) const {
  if (i > 3) return getNode(cx, cy);
  int x = 0;
  int y = 0;
  if (i == 0) {
    y = 1;
  }
  if (i == 1) {
    y = 0;
    x = 1;
  }
  if (i == 2) {
    y = -1;
  }
  if (i == 3) {
    x = -1;
  }

  return getNode(cx + x, cy + y);
}

// _____________________________________________________________________________
GridNode* GridGraph::neigh(const GridNode* n, size_t i) const {
  return neigh(n->pl().getX(), n->pl().getY(), i);
}

// _____________________________________________________________________________
void GridGraph::unSettleNd(CombNode* a) {
  openTurns(_settled[a]);
  _settled[a]->pl().setSettled(false);
  _settled.erase(a);
}

// _____________________________________________________________________________
void GridGraph::unSettleEdg(GridNode* a, GridNode* b) {
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
}

// _____________________________________________________________________________
CrossEdgPairs GridGraph::getCrossEdgPairs() const { return {}; }

// _____________________________________________________________________________
void GridGraph::addObstacle(const util::geo::Polygon<double>& obst) {
  _obstacles.push_back(obst);
  writeObstacleCost(obst);
}

// _____________________________________________________________________________
void GridGraph::writeObstacleCost(const util::geo::Polygon<double>& obst) {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto grNdA = getNode(x, y);

      for (size_t i = 0; i < maxDeg(); i++) {
        auto grNeigh = neigh(x, y, i);
        if (!grNeigh) continue;
        auto ge = getNEdg(grNdA, grNeigh);

        if (intersects(LineSegment<double>(*ge->getFrom()->pl().getGeom(),
                                           *ge->getTo()->pl().getGeom()),
                       obst) ||
            contains(LineSegment<double>(*ge->getFrom()->pl().getGeom(),
                                         *ge->getTo()->pl().getGeom()),
                     obst)) {
          ge->pl().setCost(std::numeric_limits<double>::infinity());
        }
      }
    }
  }
}

// _____________________________________________________________________________
void GridGraph::writeGeoCoursePens(const CombEdge* ce, GeoPensMap* target,
                                   double pen) {
  (*target)[ce].resize(_edgeCount);
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto grNdA = getNode(x, y);

      for (size_t i = 0; i < maxDeg(); i++) {
        auto grNeigh = neigh(x, y, i);
        if (!grNeigh) continue;
        auto ge = getNEdg(grNdA, grNeigh);

        double d = std::numeric_limits<double>::infinity();

        for (auto orE : ce->pl().getChilds()) {
          double dLoc =
              dist(*orE->pl().getGeom(),
                   LineSegment<double>(*ge->getFrom()->pl().getGeom(),
                                       *ge->getTo()->pl().getGeom())) /
              getCellSize();
          if (dLoc < d) d = dLoc;
        }

        d *= pen * d;

        (*target)[ce][ge->pl().getId()] = d;
      }
    }
  }
}

// _____________________________________________________________________________
void GridGraph::settleEdg(GridNode* a, GridNode* b, CombEdge* e) {
  if (a == b) return;

  // this closes the grid edge
  auto ge = getNEdg(a, b);

  assert(ge);

  addResEdg(ge, e);

  // this closes both nodes
  // a close means that all major edges reaching this node are closed
  closeTurns(a);
  closeTurns(b);
}

// _____________________________________________________________________________
void GridGraph::addResEdg(GridEdge* ge, CombEdge* ce) {
  assert(_resEdgs.count(ge) == 0 || _resEdgs.find(ge)->second == 0);
  ge->pl().addResEdge();
  _resEdgs[ge] = ce;
}

// _____________________________________________________________________________
CombEdge* GridGraph::getResEdg(GridEdge* ge) {
  if (!ge) return 0;
  if (_resEdgs.count(ge)) return _resEdgs.find(ge)->second;
  return 0;
}

// _____________________________________________________________________________
GridEdge* GridGraph::getNEdg(const GridNode* a, const GridNode* b) const {
  if (!a || !b) return 0;

  int aa = 1 + (int)a->pl().getX() - (int)b->pl().getX();
  int bb = 1 + (int)a->pl().getY() - (int)b->pl().getY();

  size_t dir = 0;
  size_t d = aa * 3 + bb;

  if (d == 1)
    dir = 1;
  else if (d == 3)
    dir = 0;
  else if (d == 5)
    dir = 2;
  else if (d == 7)
    dir = 3;
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
void GridGraph::getSettledAdjEdgs(GridNode* n, CombEdge* outgoing[8]) {
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t i = 0; i < maxDeg(); i++) {
    outgoing[i] = 0;
    auto p = n->pl().getPort(i);
    auto neighbor = neigh(x, y, i);

    if (!neighbor) continue;

    auto neighP = neighbor->pl().getPort((i + maxDeg() / 2) % maxDeg());
    auto e = getEdg(p, neighP);
    auto f = getEdg(neighP, p);
    auto resEdg = getResEdg(e);
    if (!resEdg) resEdg = getResEdg(f);

    if (resEdg) outgoing[i] = resEdg;
  }
}

// _____________________________________________________________________________
const Penalties& GridGraph::getPens() const { return _c; }

// _____________________________________________________________________________
NodeCost GridGraph::nodeBendPen(GridNode* n, CombEdge* e) {
  NodeCost addC;

  CombEdge* out[maxDeg()];
  getSettledAdjEdgs(n, out);

  for (size_t i = 0; i < maxDeg(); i++) {
    if (!out[i]) continue;
    for (auto lo : e->pl().getChilds().front()->pl().getLines()) {
      // TODO: turn restrictions, if there is actually no connection
      // between the lines on the edges, dont penalize!!
      if (out[i]->pl().getChilds().front()->pl().hasLine(lo.line)) {
        for (size_t j = 0; j < maxDeg(); j++) {
          addC[j] += getBendPen(i, j);
        }

        // only count a single matching line
        break;
      }
    }
  }

  return addC;
}

// _____________________________________________________________________________
double GridGraph::getBendPen(size_t i, size_t j) const {
  return _bendCosts[ang(i, j)];
}

// _____________________________________________________________________________
size_t GridGraph::ang(size_t i, size_t j) const {
  // determine angle between port i and j
  int ang = (4 + (i - j)) % 4;
  if (ang > 2) ang = 4 - ang;
  ang = ang % 2;

  int d = (int)(i) - (int)(j);
  int deg = abs((((d + 2) % 4) + 4) % 4 - 2) % 2;

  assert(ang == deg);

  return ang;
}

// _____________________________________________________________________________
NodeCost GridGraph::spacingPen(GridNode* nd, CombNode* origNd, CombEdge* edg) {
  NodeCost addC;

  CombEdge* out[maxDeg()];
  getSettledAdjEdgs(nd, out);

  for (size_t i = 0; i < maxDeg(); i++) {
    if (!out[i]) continue;

    // this is the number of edges that will occur between the currently checked
    // edge and the inserted edge, in clockwise and counter-clockwise dir
    int32_t dCw = origNd->pl().getEdgeOrdering().dist(out[i], edg);
    int32_t dCCw = origNd->pl().getEdgeOrdering().dist(edg, out[i]);

    for (int j = 1; j < dCw; j++) {
      addC[(i + j) % maxDeg()] = -1.0 * std::numeric_limits<double>::max();
    }

    for (int j = 1; j < dCCw; j++) {
      addC[(i + (maxDeg() - j)) % maxDeg()] =
          -1.0 * std::numeric_limits<double>::max();
    }
  }

  return addC;
}

// _____________________________________________________________________________
NodeCost GridGraph::topoBlockPen(GridNode* nd, CombNode* origNd,
                                 CombEdge* edg) {
  CombEdge* outgoing[maxDeg()];
  NodeCost addC;
  getSettledAdjEdgs(nd, outgoing);

  // topological blocking
  for (size_t i = 0; i < maxDeg(); i++) {
    if (!outgoing[i]) continue;

    for (size_t j = i + 1; j < i + maxDeg(); j++) {
      if (!outgoing[j % maxDeg()]) continue;
      if (outgoing[j % maxDeg()] == outgoing[i]) break;

      int da = origNd->pl().getEdgeOrdering().dist(outgoing[i], edg);
      int db = origNd->pl().getEdgeOrdering().dist(outgoing[j % maxDeg()], edg);

      if (db < da) {
        // edge does not lie in this segment, block it!
        for (size_t x = i + 1; x < j; x++) {
          addC[x % maxDeg()] = -1.0 * std::numeric_limits<double>::max();
        }
      }
    }
  }
  return addC;
}

// _____________________________________________________________________________
void GridGraph::addCostVec(GridNode* n, const NodeCost& addC) {
  for (size_t i = 0; i < maxDeg(); i++) {
    auto p = n->pl().getPort(i);

    if (addC[i] < -1) {
      getEdg(p, n)->pl().close();
      getEdg(n, p)->pl().close();
    } else {
      getEdg(p, n)->pl().setCost(getEdg(p, n)->pl().rawCost() + addC[i]);
      getEdg(n, p)->pl().setCost(getEdg(n, p)->pl().rawCost() + addC[i]);
    }
  }
}

// _____________________________________________________________________________
void GridGraph::writeInitialCosts() {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto n = getNode(x, y);
      for (size_t i = 0; i < maxDeg(); i++) {
        auto port = n->pl().getPort(i);
        auto neighbor = neigh(x, y, i);

        if (!neighbor || !port) continue;

        auto oPort = neighbor->pl().getPort((i + maxDeg() / 2) % maxDeg());
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
std::priority_queue<Candidate> GridGraph::getGridNdCands(const DPoint& p,
                                                         double maxD) const {
  std::priority_queue<Candidate> ret;
  std::set<GridNode*> neigh;
  DBox b(DPoint(p.getX() - maxD, p.getY() - maxD),
         DPoint(p.getX() + maxD, p.getY() + maxD));
  _grid.get(b, &neigh);

  for (auto n : neigh) {
    if (n->pl().isClosed() || n->pl().isSettled()) continue;
    double d = dist(*n->pl().getGeom(), p);
    if (d < maxD) ret.push(Candidate(n, d));
  }

  return ret;
}

// _____________________________________________________________________________
const Grid<GridNode*, Point, double>& GridGraph::getGrid() const {
  return _grid;
}

// _____________________________________________________________________________
double GridGraph::heurCost(int64_t xa, int64_t ya, int64_t xb,
                           int64_t yb) const {
  int dx = labs(xb - xa);
  int dy = labs(yb - ya);

  // Alternative: use chebyshev distance heuristic
  // double minHops = std::max(dx, dy);

  // double heurECost =
  // (std::min(_c.verticalPen, std::min(_c.horizontalPen, _c.diagonalPen)));

  // return minHops * (heurECost + _heurHopCost) - _heurHopCost;

  double edgCost = ((_c.horizontalPen + _heurHopCost) * dx +
                    (_c.verticalPen + _heurHopCost) * dy);

  // we have to do at least one turn, which can only be a 90 degree turn
  if (dx != 0 && dy != 0) edgCost += _c.p_90;

  // we always count one heurHopCost too much, subtract it at the end!
  return edgCost - _heurHopCost;
}

// _____________________________________________________________________________
const util::graph::Dijkstra::HeurFunc<GridNodePL, GridEdgePL, float>*
GridGraph::getHeur(const std::set<GridNode*>& to) const {
  return new GridGraphHeur(this, to);
}

// _____________________________________________________________________________
void GridGraph::openTurns(GridNode* n) {
  if (!n->pl().isClosed()) return;

  // open all non-sink inner edges
  for (size_t i = 0; i < maxDeg(); i++) {
    auto portA = n->pl().getPort(i);
    if (portA->getDeg() == maxDeg()) continue;
    for (size_t j = i + 1; j < maxDeg(); j++) {
      auto portB = n->pl().getPort(j);
      if (portB->getDeg() == maxDeg()) continue;
      auto e = getEdg(portA, portB);
      auto f = getEdg(portB, portA);

      e->pl().open();
      f->pl().open();
    }
  }

  n->pl().setClosed(false);
}

// _____________________________________________________________________________
void GridGraph::closeTurns(GridNode* n) {
  if (n->pl().isClosed()) return;

  // close all non-sink inner edges
  for (size_t i = 0; i < maxDeg(); i++) {
    auto portA = n->pl().getPort(i);
    for (size_t j = i + 1; j < maxDeg(); j++) {
      auto portB = n->pl().getPort(j);
      auto e = getEdg(portA, portB);
      auto f = getEdg(portB, portA);

      e->pl().close();
      f->pl().close();
    }
  }

  n->pl().setClosed(true);
}

// _____________________________________________________________________________
void GridGraph::openSinkTo(GridNode* n, double cost) {
  for (size_t i = 0; i < maxDeg(); i++) {
    getEdg(n->pl().getPort(i), n)->pl().open();
    getEdg(n->pl().getPort(i), n)->pl().setCost(cost);
  }
}

// _____________________________________________________________________________
void GridGraph::closeSinkTo(GridNode* n) {
  for (size_t i = 0; i < maxDeg(); i++) {
    getEdg(n->pl().getPort(i), n)->pl().close();
    getEdg(n->pl().getPort(i), n)->pl().setCost(INF);
  }
}

// _____________________________________________________________________________
void GridGraph::openSinkFr(GridNode* n, double cost) {
  for (size_t i = 0; i < maxDeg(); i++) {
    getEdg(n, n->pl().getPort(i))->pl().open();
    getEdg(n, n->pl().getPort(i))->pl().setCost(cost);
  }
}

// _____________________________________________________________________________
void GridGraph::closeSinkFr(GridNode* n) {
  for (size_t i = 0; i < maxDeg(); i++) {
    getEdg(n, n->pl().getPort(i))->pl().close();
    getEdg(n, n->pl().getPort(i))->pl().setCost(INF);
  }
}

// _____________________________________________________________________________
GridNode* GridGraph::getSettled(const CombNode* cnd) const {
  if (_settled.count(cnd)) return _settled.find(cnd)->second;
  return 0;
}

// _____________________________________________________________________________
std::set<GridNode*> GridGraph::getGrNdCands(CombNode* n, double maxDis) {
  std::set<GridNode*> tos;
  if (!isSettled(n)) {
    auto cands = getGridNdCands(*n->pl().getGeom(), maxDis);

    while (!cands.empty()) {
      if (!cands.top().n->pl().isClosed()) tos.insert(cands.top().n);
      cands.pop();
    }
  } else {
    tos.insert(_settled.find(n)->second);
  }

  return tos;
}

// _____________________________________________________________________________
void GridGraph::settleNd(GridNode* n, CombNode* cn) {
  _settled[cn] = n;
  n->pl().setSettled(true);
}

// _____________________________________________________________________________
bool GridGraph::isSettled(const CombNode* cn) {
  return _settled.find(cn) != _settled.end();
}

// _____________________________________________________________________________
GridNode* GridGraph::getGrNdById(size_t id) const { return _nds[id]; }

// _____________________________________________________________________________
const GridEdge* GridGraph::getGrEdgById(std::pair<size_t, size_t> id) const {
  return getEdg(_nds[id.first], _nds[id.second]);
}

// _____________________________________________________________________________
GridNode* GridGraph::writeNd(size_t x, size_t y) {
  double xPos = _bbox.getLowerLeft().getX() + x * _cellSize;
  double yPos = _bbox.getLowerLeft().getY() + y * _cellSize;

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
      double pen = getBendPen(i, j);

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
double GridGraph::ndMovePen(const CombNode* cbNd, const GridNode* grNd) const {
  // the move penalty has to be at least the max cost of saving a single
  // grid hop - otherwise we could move the node closer and closer to the
  // other node without ever increasing the total cost.

  // additional penalty per grid move
  // TODO: make configurable
  double PEN = 0.5;

  double c_0 = _c.p_45 - _c.p_135;
  double penPerGrid = PEN + c_0 + fmax(_c.diagonalPen, _c.horizontalPen);

  // real distance between grid node and input node
  double d = dist(*cbNd->pl().getGeom(), *grNd->pl().getGeom());

  // distance normalized to grid length
  double gridD = d / getCellSize();

  // and multiplied per grid hop penalty
  return gridD * penPerGrid;
}

// _____________________________________________________________________________
void GridGraph::reset() {
  _settled.clear();
  _resEdgs.clear();
  for (auto n : *getNds()) {
    for (auto e : n->getAdjListOut()) e->pl().reset();
    if (!n->pl().isSink()) continue;
    openTurns(n);
    closeSinkFr(n);
    closeSinkTo(n);
  }

  writeInitialCosts();
  reWriteObstCosts();
}

// _____________________________________________________________________________
void GridGraph::reWriteObstCosts() {
  for (const auto& obst : _obstacles) {
    writeObstacleCost(obst);
  }
}

// _____________________________________________________________________________
double GridGraph::getCellSize() const { return _cellSize; }
