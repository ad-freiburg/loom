// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include "octi/gridgraph/GridGraph.h"
#include "octi/gridgraph/NodeCost.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Node.h"

using namespace octi::gridgraph;
using octi::gridgraph::NodeCost;
using octi::gridgraph::GridGraph;
using util::geo::DBox;
using util::geo::DPoint;
using util::geo::dist;

double INF = std::numeric_limits<float>::infinity();

// _____________________________________________________________________________
GridGraph::GridGraph(const DBox& bbox, double cellSize, double spacer,
                     const Penalties& pens)
    : _bbox(bbox),
      _c(pens),
      _grid(cellSize, cellSize, bbox, false),
      _cellSize(cellSize),
      _spacer(spacer) {
  assert(_c.p_0 <= _c.p_135);
  assert(_c.p_135 <= _c.p_90);
  assert(_c.p_90 <= _c.p_45);

  double c_0 = _c.p_45 - _c.p_135;
  double c_135 = _c.p_45;
  double c_90 = _c.p_45 - _c.p_135 + _c.p_90;
  double c_45 = c_0 + c_135;

  _heurECost =
      (std::min(_c.verticalPen, std::min(_c.horizontalPen, _c.diagonalPen)));
  _heurHopCost = _c.p_45 - _c.p_135;

  // std::cerr << "C0: " << c_0 << std::endl;
  // std::cerr << "C135: " << c_135 << std::endl;
  // std::cerr << "C90: " << c_90 << std::endl;
  // std::cerr << "C45: " << c_45 << std::endl;

  // cut off illegal spacer values
  if (spacer > cellSize / 2) spacer = cellSize / 2;

  // write nodes
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      writeNd(x, y);
    }
  }

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      GridNode* center = getNode(x, y);

      for (size_t p = 0; p < 8; p++) {
        GridNode* from = center->pl().getPort(p);
        GridNode* toN = getNeighbor(x, y, p);
        if (from != 0 && toN != 0) {
          GridNode* to = toN->pl().getPort((p + 4) % 8);
          auto e = new GridEdge(from, to, GridEdgePL(9, false));
          from->addEdge(e);
          to->addEdge(e);
        }
      }
    }
  }

  writeInitialCosts();
}

// _____________________________________________________________________________
GridNode* GridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  return _nds[_grid.getYHeight() * 9 * x + y * 9];
}

// _____________________________________________________________________________
GridNode* GridGraph::getNeighbor(size_t cx, size_t cy, size_t i) const {
  int8_t x = 1;
  if (i % 4 == 0) x = 0;
  if (i > 4) x = -1;

  int8_t y = 1;
  if (i == 2 || i == 6) y = 0;
  if (i == 3 || i == 4 || i == 5) y = -1;

  return getNode(cx + x, cy + y);
}

// _____________________________________________________________________________
void GridGraph::unSettleNd(CombNode* a) {
  openNodeTurns(_settled[a]);
  _settled[a]->pl().setSettled(false);
  _settled.erase(a);
}

// _____________________________________________________________________________
void GridGraph::unSettleEdg(GridNode* a, GridNode* b) {
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
    openNodeTurns(a);
  }
  if (!b->pl().isSettled()) {
    openNodeTurns(b);
  }

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
void GridGraph::writeGeoCoursePens(const CombEdge* ce) {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto grNdA = getNode(x, y);

      for (size_t i = 0; i < 8; i++) {
        auto neigh = getNeighbor(x, y, i);
        if (!neigh) continue;
        auto ge = getNEdg(grNdA, neigh);

        double d = std::numeric_limits<double>::infinity();

        for (auto orE : ce->pl().getChilds()) {
          double dLoc = util::geo::dist(*orE->pl().getGeom(),
                                        util::geo::LineSegment<double>(
                                            *ge->getFrom()->pl().getGeom(),
                                            *ge->getTo()->pl().getGeom())) /
                        getCellSize();
          if (dLoc < d) d = dLoc;
        }

        d *= 2 * d;

        ge->pl().setGeoCourseCost(d);
      }
    }
  }
}

// _____________________________________________________________________________
void GridGraph::clearGeoCoursePens() {
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto grNdA = getNode(x, y);

      for (size_t i = 0; i < 8; i++) {
        auto neigh = getNeighbor(x, y, i);
        if (!neigh) continue;
        auto ge = getNEdg(grNdA, neigh);

        ge->pl().clearGeoCourseCost();
      }
    }
  }
}

// _____________________________________________________________________________
void GridGraph::settleEdg(GridNode* a, GridNode* b, CombEdge* e) {
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
  // auto gf = getNEdg(b, a);

  addResEdg(ge, e);
  // addResEdg(gf, e);

  // this closes both nodes
  // a close means that all major edges reaching this node are closed
  closeNodeTurns(a);
  closeNodeTurns(b);

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
  if (!a) return 0;
  if (!b) return 0;

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

  if (a->pl().getPort(dir) && b->pl().getPort((dir + 4) % 8)) {
    return const_cast<GridEdge*>(
        getEdg(a->pl().getPort(dir), b->pl().getPort((dir + 4) % 8)));
  }

  return 0;
}

// _____________________________________________________________________________
void GridGraph::getSettledAdjEdgs(GridNode* n, CombEdge* outgoing[8]) {
  size_t x = n->pl().getX();
  size_t y = n->pl().getY();

  for (size_t i = 0; i < 8; i++) {
    outgoing[i] = 0;
    auto p = n->pl().getPort(i);
    auto neigh = getNeighbor(x, y, i);

    if (!neigh) continue;

    auto neighP = neigh->pl().getPort((i + 4) % 8);
    auto e = getEdg(p, neighP);
    auto f = getEdg(neighP, p);
    auto resEdg = getResEdg(e);
    if (!resEdg) resEdg = getResEdg(f);

    if (resEdg) outgoing[i] = resEdg;
  }
}

// _____________________________________________________________________________
const Penalties& GridGraph::getPenalties() const { return _c; }

// _____________________________________________________________________________
NodeCost GridGraph::nodeBendPenalty(GridNode* n, CombEdge* e) {
  NodeCost addC;

  // TODO: same code as in write nd
  double c_0 = _c.p_45 - _c.p_135;
  double c_135 = _c.p_45;
  double c_90 = _c.p_45 - _c.p_135 + _c.p_90;
  double c_45 = c_0 + c_135;

  CombEdge* out[8];
  getSettledAdjEdgs(n, out);

  for (int i = 0; i < 8; i++) {
    if (!out[i]) continue;
    for (auto ro : e->pl().getChilds().front()->pl().getRoutes()) {
      // TODO: turn restrictions, if there is actually no connection
      // between the lines on the edges, dont penalize!!
      if (out[i]->pl().getChilds().front()->pl().hasRoute(ro.route)) {
        for (int j = 0; j < 8; j++) {
          // determine angle between port i and j

          int ang = (8 + (i - j)) % 8;
          if (ang > 4) ang = 8 - ang;
          ang = ang % 4;

          double mult = 1;

          // write corresponding cost to addC[j]
          if (ang == 0) addC[j] += mult * c_0;
          if (ang == 1) addC[j] += mult * c_45;
          if (ang == 2) addC[j] += mult * c_90;
          if (ang == 3) addC[j] += mult * c_135;
        }

        // only count a single matching line
        break;
      }
    }
  }

  return addC;
}

// _____________________________________________________________________________
NodeCost GridGraph::spacingPenalty(GridNode* nd, CombNode* origNd,
                                   CombEdge* edg) {
  NodeCost addC;

  CombEdge* out[8];
  getSettledAdjEdgs(nd, out);

  for (size_t i = 0; i < 8; i++) {
    if (!out[i]) continue;

    // this is the number of edges that will occur between the currently checked
    // edge and the inserted edge, in clockwise and counter-clockwise dir
    int32_t dCw = origNd->pl().getEdgeOrdering().dist(out[i], edg);
    int32_t dCCw = origNd->pl().getEdgeOrdering().dist(edg, out[i]);

    for (int j = 1; j < dCw; j++) {
      addC[(i + j) % 8] = -1.0 * std::numeric_limits<double>::max();
    }

    for (int j = 1; j < dCCw; j++) {
      addC[(i + (8 - j)) % 8] = -1.0 * std::numeric_limits<double>::max();
    }
  }

  return addC;
}

// _____________________________________________________________________________
NodeCost GridGraph::topoBlockPenalty(GridNode* nd, CombNode* origNd,
                                     CombEdge* edg) {
  CombEdge* outgoing[8];
  NodeCost addC;
  getSettledAdjEdgs(nd, outgoing);

  // topological blocking
  for (size_t i = 0; i < 8; i++) {
    if (!outgoing[i]) continue;

    for (size_t j = i + 1; j < i + 8; j++) {
      if (!outgoing[j % 8]) continue;
      if (outgoing[j % 8] == outgoing[i]) break;

      int da = origNd->pl().getEdgeOrdering().dist(outgoing[i], edg);
      int db = origNd->pl().getEdgeOrdering().dist(outgoing[j % 8], edg);

      if (db < da) {
        // edge does not lie in this segment, block it!
        for (size_t x = i + 1; x < j; x++) {
          addC[x % 8] = -1.0 * std::numeric_limits<double>::max();
        }
      }
    }
  }
  return addC;
}

// _____________________________________________________________________________
void GridGraph::addCostVector(GridNode* n, const NodeCost& addC) {
  for (size_t i = 0; i < 8; i++) {
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
      for (size_t i = 0; i < 8; i++) {
        auto port = n->pl().getPort(i);
        auto neigh = getNeighbor(x, y, i);

        if (!neigh || !port) continue;

        auto oPort = neigh->pl().getPort((i + 4) % 8);
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
  double minHops = std::max(labs(xb - xa), labs(yb - ya));
  return minHops * (_heurECost + _heurHopCost) - _heurHopCost;
}

// _____________________________________________________________________________
void GridGraph::openNodeTurns(GridNode* n) {
  if (!n->pl().isClosed()) return;

  // open all non-sink inner edges
  for (size_t i = 0; i < 8; i++) {
    auto portA = n->pl().getPort(i);
    if (portA->getDeg() == 8) continue;
    for (size_t j = i + 1; j < 8; j++) {
      auto portB = n->pl().getPort(j);
      if (portB->getDeg() == 8) continue;
      auto e = getEdg(portA, portB);
      auto f = getEdg(portB, portA);

      e->pl().open();
      f->pl().open();
    }
  }

  n->pl().setClosed(false);
}

// _____________________________________________________________________________
void GridGraph::closeNodeTurns(GridNode* n) {
  if (n->pl().isClosed()) return;

  // close all non-sink inner edges
  for (size_t i = 0; i < 8; i++) {
    auto portA = n->pl().getPort(i);
    for (size_t j = i + 1; j < 8; j++) {
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
void GridGraph::openNodeSinkTo(GridNode* n, double cost) {
  for (size_t i = 0; i < 8; i++) {
    getEdg(n->pl().getPort(i), n)->pl().open();
    getEdg(n->pl().getPort(i), n)->pl().setCost(cost);
  }
}

// _____________________________________________________________________________
void GridGraph::closeNodeSinkTo(GridNode* n) {
  for (size_t i = 0; i < 8; i++) {
    getEdg(n->pl().getPort(i), n)->pl().close();
    getEdg(n->pl().getPort(i), n)->pl().setCost(INF);
  }
}

// _____________________________________________________________________________
void GridGraph::openNodeSinkFr(GridNode* n, double cost) {
  for (size_t i = 0; i < 8; i++) {
    getEdg(n, n->pl().getPort(i))->pl().open();
    getEdg(n, n->pl().getPort(i))->pl().setCost(cost);
  }
}

// _____________________________________________________________________________
void GridGraph::closeNodeSinkFr(GridNode* n) {
  for (size_t i = 0; i < 8; i++) {
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
    n->addEdge(e);
    nn->addEdge(e);

    e = new GridEdge(nn, n, GridEdgePL(INF, true, false));
    n->addEdge(e);
    nn->addEdge(e);
  }

  // in-node connections
  for (size_t i = 0; i < 8; i++) {
    for (size_t j = i + 1; j < 8; j++) {
      int d = (int)(i) - (int)(j);
      size_t deg = abs((((d + 4) % 8) + 8) % 8 - 4);
      double pen = c_0;

      if (deg == 0) pen = c_0;
      if (deg == 1) pen = c_45;
      if (deg == 2) pen = c_90;
      if (deg == 3) pen = c_135;

      if (x == 0 && (i == 5 || i == 6 || i == 7)) pen = INF;
      if (y == 0 && (i == 0 || i == 7 || i == 1)) pen = INF;
      if (x == _grid.getXWidth() - 1 && (i == 1 || i == 2 || i == 3)) pen = INF;
      if (y == _grid.getYHeight() - 1 && (i == 3 || i == 4 || i == 5))
        pen = INF;

      auto e = new GridEdge(n->pl().getPort(i), n->pl().getPort(j),
                            GridEdgePL(pen, true));
      e->getFrom()->addEdge(e);
      e->getTo()->addEdge(e);

      e = new GridEdge(n->pl().getPort(j), n->pl().getPort(i),
                       GridEdgePL(pen, true));
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
    openNodeTurns(n);
    closeNodeSinkFr(n);
    closeNodeSinkTo(n);
  }

  writeInitialCosts();
}

// _____________________________________________________________________________
double GridGraph::getCellSize() const { return _cellSize; }
