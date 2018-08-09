// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include "octi/gridgraph/GridGraph.h"
#include "util/graph/Node.h"

using namespace octi::gridgraph;

double INF = std::numeric_limits<double>::infinity();

// _____________________________________________________________________________
GridGraph::GridGraph(const util::geo::DBox& bbox, double cellSize,
                     const Penalties& pens)
    : _bbox(bbox), _c(pens), _grid(cellSize, cellSize, bbox) {
  assert(_c.p_0 < _c.p_135);
  assert(_c.p_135 < _c.p_90);
  assert(_c.p_90 < _c.p_45);

  double c_0 = _c.p_45 - _c.p_135;
  double c_135 = _c.p_45;
  double c_90 = _c.p_45 - _c.p_135 + _c.p_90;

  double spacer = cellSize / 10.0;

  // write nodes
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      double xPos = bbox.getLowerLeft().getX() + x * cellSize;
      double yPos = bbox.getLowerLeft().getY() + y * cellSize;
      GridNode* n = addNd(DPoint(xPos, yPos));
      _grid.add(x, y, n);
      n->pl().setXY(x, y);
      n->pl().setParent(n);

      GridNode* n1 = addNd(DPoint(xPos - spacer, yPos));
      n1->pl().setParent(n);
      n->pl().setPort(6, n1);
      addEdg(n, n1, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n1->pl().getGeom()),
                           INF, true, false));
      addEdg(n1, n, GridEdgePL(util::geo::PolyLine<double>(*n1->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));

      GridNode* n2 =
          addNd(DPoint(xPos - spacer, yPos - spacer));
      n2->pl().setParent(n);
      n->pl().setPort(5, n2);
      addEdg(n, n2, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n2->pl().getGeom()),
                           INF, true, false));
      addEdg(n2, n, GridEdgePL(util::geo::PolyLine<double>(*n2->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));

      GridNode* n3 = addNd(DPoint(xPos, yPos - spacer));
      n3->pl().setParent(n);
      n->pl().setPort(4, n3);
      addEdg(n, n3, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n3->pl().getGeom()),
                           INF, true, false));
      addEdg(n3, n, GridEdgePL(util::geo::PolyLine<double>(*n3->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));

      GridNode* n4 = addNd(DPoint(xPos + spacer, yPos));
      n4->pl().setParent(n);
      n->pl().setPort(2, n4);
      addEdg(n, n4, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n4->pl().getGeom()),
                           INF, true, false));
      addEdg(n4, n, GridEdgePL(util::geo::PolyLine<double>(*n4->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));

      GridNode* n5 =
          addNd(DPoint(xPos + spacer, yPos + spacer));
      n5->pl().setParent(n);
      n->pl().setPort(1, n5);
      addEdg(n, n5, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n5->pl().getGeom()),
                           INF, true, false));
      addEdg(n5, n, GridEdgePL(util::geo::PolyLine<double>(*n5->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));

      GridNode* n6 = addNd(DPoint(xPos, yPos + spacer));
      n6->pl().setParent(n);
      n->pl().setPort(0, n6);
      addEdg(n, n6, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n6->pl().getGeom()),
                           INF, true, false));
      addEdg(n6, n, GridEdgePL(util::geo::PolyLine<double>(*n6->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));

      GridNode* n7 =
          addNd(DPoint(xPos - spacer, yPos + spacer));
      n7->pl().setParent(n);
      n->pl().setPort(7, n7);
      addEdg(n, n7, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n7->pl().getGeom()),
                           INF, true, false));
      addEdg(n7, n, GridEdgePL(util::geo::PolyLine<double>(*n7->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));
      GridNode* n8 =
          addNd(DPoint(xPos + spacer, yPos - spacer));
      n8->pl().setParent(n);
      n->pl().setPort(3, n8);
      addEdg(n, n8, GridEdgePL(util::geo::PolyLine<double>(*n->pl().getGeom(),
                                                       *n8->pl().getGeom()),
                           INF, true, false));
      addEdg(n8, n, GridEdgePL(util::geo::PolyLine<double>(*n8->pl().getGeom(),
                                                       *n->pl().getGeom()),
                           INF, true, false));

      // in-node connections
      for (size_t i = 0; i < 8; i++) {
        for (size_t j = i + 1; j < 8; j++) {
          int d = (int)(i) - (int)(j);
          size_t deg = abs((((d + 4) % 8) + 8) % 8 - 4);
          double pen = c_0;

          if (deg == 1) continue;
          if (deg == 2) pen = c_90;
          if (deg == 3) pen = c_135;
          addEdg(n->pl().getPort(i), n->pl().getPort(j),
                 GridEdgePL(util::geo::PolyLine<double>(
                            *n->pl().getPort(i)->pl().getGeom(),
                            *n->pl().getPort(j)->pl().getGeom()),
                        pen, true));
          addEdg(n->pl().getPort(j), n->pl().getPort(i),
                 GridEdgePL(util::geo::PolyLine<double>(
                            *n->pl().getPort(j)->pl().getGeom(),
                            *n->pl().getPort(i)->pl().getGeom()),
                        pen, true));
        }
      }
    }
  }

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      GridNode* center = getNode(x, y);
      if (center == 0) continue;

      for (size_t p = 0; p < 8; p++) {
        GridNode* from = center->pl().getPort(p);
        GridNode* toN = getNeighbor(x, y, p);
        if (from != 0 && toN != 0) {
          GridNode* to = toN->pl().getPort((p + 4) % 8);
          addEdg(from, to,
                 GridEdgePL(util::geo::PolyLine<double>(*from->pl().getGeom(),
                                                    *to->pl().getGeom()),
                        0, false));
        }
      }
    }
  }

  writeInitialCosts();
}

// _____________________________________________________________________________
GridNode* GridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  std::set<GridNode*> r;
  _grid.get(x, y, &r);

  if (r.size()) return *r.begin();
  return 0;
}

// _____________________________________________________________________________
std::pair<size_t, size_t> GridGraph::getNodeCoords(GridNode* n) const {
  return std::pair<size_t, size_t>(n->pl().getX(), n->pl().getY());
}

// _____________________________________________________________________________
GridNode* GridGraph::getNeighbor(size_t cx, size_t cy,
                                             size_t i) const {
  int8_t x = 1;
  if (i % 4 == 0) x = 0;
  if (i > 4) x = -1;

  int8_t y = 1;
  if (i == 2 || i == 6) y = 0;
  if (i == 3 || i == 4 || i == 5) y = -1;

  return getNode(cx + x, cy + y);
}

// _____________________________________________________________________________
void GridGraph::balanceEdge(GridNode* a, GridNode* b) {
  if (a == b) return;
  size_t dir = 0;
  for (; dir < 8; dir++) {
    if (getEdg(a->pl().getPort(dir), b->pl().getPort((dir + 4) % 8))) {
      break;
    }
  }

  auto xy = getNodeCoords(a);
  size_t x = xy.first;
  size_t y = xy.second;

  getNEdge(a, b)->pl().setCost(INF);
  getNEdge(b, a)->pl().setCost(INF);

  closeNode(a);
  closeNode(b);

  if (dir == 1 || dir == 3 || dir == 5 || dir == 7) {
    auto na = getNeighbor(x, y, (dir + 7) % 8);
    auto nb = getNeighbor(x, y, (dir + 1) % 8);

    if (na && nb) {
      auto e = getNEdge(na, nb);
      auto f = getNEdge(nb, na);

      e->pl().setCost(INF);
      f->pl().setCost(INF);
    }
  }
}

// _____________________________________________________________________________
GridEdge* GridGraph::getOtherEdge(GridEdge* e) {
  return getEdg(e->getTo(), e->getFrom());
}

// _____________________________________________________________________________
GridEdge* GridGraph::getNEdge(GridNode* a, GridNode* b) {
  for (size_t i = 0; i < 8; i++) {
    if (a && b && a->pl().getPort(i) && b->pl().getPort((i + 4) % 8)) {
      auto e = getEdg(a->pl().getPort(i), b->pl().getPort((i + 4) % 8));
      if (e) return e;
    }
  }

  return 0;
}

// _____________________________________________________________________________
void GridGraph::getSettledOutgoingEdges(GridNode* n, CombEdge* outgoing[8]) {
  auto xy = getNodeCoords(n);
  size_t x = xy.first;
  size_t y = xy.second;

  // if some outgoing edge is taken, dont put new edge next to it
  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    auto neigh = getNeighbor(x, y, i);

    if (neigh &&
        getEdg(port, neigh->pl().getPort((i + 4) % 8))
                ->pl()
                .getResEdges()
                .size() > 0) {
      outgoing[i] = *getEdg(port, neigh->pl().getPort((i + 4) % 8))
                         ->pl()
                         .getResEdges()
                         .begin();
    } else {
      outgoing[i] = 0;
    }
  }
}

// _____________________________________________________________________________
void GridGraph::spacingPenalty(GridNode* n, CombNode* origNode, CombEdge* e,
                               double* addC) {
  auto xy = getNodeCoords(n);
  size_t x = xy.first;
  size_t y = xy.second;
  assert(getNode(x, y) == n);

  std::cerr << "Checking " << origNode->pl().toString() << std::endl;

  int origEdgeNumber = origNode->getAdjList().size();
  size_t optimDistance = (8 / origEdgeNumber) - 1;

  if (!origNode->pl().getEdgeOrdering().has(e)) {
    std::cerr << "Warning: tried to balance edge " << e << " in node "
              << origNode << ", but the edge does not appear there."
              << std::endl;
    return;
  }

  std::cerr << std::endl;
  std::cerr << "Orig edge number = " << origEdgeNumber << ", optim distance is "
            << optimDistance << std::endl;

  CombEdge* outgoing[8];
  getSettledOutgoingEdges(n, outgoing);

  std::cerr << "Edge distribution: ";
  if (origNode->pl().getParent()->pl().getStops().size())
    std::cerr << "in "
              << origNode->pl().getParent()->pl().getStops().front().name;
  std::cerr << " ";
  for (size_t i = 0; i < 8; i++) {
    std::cerr << outgoing[i] << ", ";
  }
  std::cerr << std::endl;

  std::cerr << "Ordered comb edges in orig node: ";
  std::cerr << origNode->pl().getEdgeOrdering().toString(origNode);
  std::cerr << std::endl;

  for (size_t i = 0; i < 8; i++) {
    if (!outgoing[i]) continue;

    // this is the number of edges that will occur between the currently checked
    // edge and the inserted edge, in clockwise and counter-clockwise dir
    int32_t dCw = origNode->pl().getEdgeOrdering().dist(outgoing[i], e) - 1;
    int32_t dCCw = origNode->pl().getEdgeOrdering().dist(e, outgoing[i]) - 1;

    // dd and ddd are the optimal distances between outgoing[i] and e, based on
    // the total number
    // of edges in this node
    int dd = ((((dCw + 1) + dCw) % 8) * optimDistance) % 8;
    int ddd = (6 - dd) % 8;

    std::cerr << "Distance between the inserted edge (" << e << ") and edge at "
              << i << " (" << outgoing[i] << ") is " << dCw << " (cw) and "
              << dCCw << " (ccw)"
              << ", optim distance between them is +" << dd << " and -" << ddd
              << std::endl;

    double pen = _c.p_45 * 2 - 1;

    for (int j = 1; dd != 0 && j <= dd + 1; j++) {
      if (addC[(i + j) % 8] < -1) continue;
      addC[(i + j) % 8] += pen * (1.0 - (j - 1.0) / (dd));
    }

    for (int j = 1; ddd != 0 && j <= ddd + 1; j++) {
      if (addC[(i + (8 - j)) % 8] < -1) continue;
      addC[(i + (8 - j)) % 8] += pen * (1.0 - (j - 1.0) / (ddd));
    }

    // negative cost here means that the edge is going to be closed
    addC[i] = -1.0 * std::numeric_limits<double>::max();

    for (int j = 1; j <= dCw; j++) {
      addC[(i + j) % 8] = -1.0 * std::numeric_limits<double>::max();
    }

    for (int j = 1; j <= dCCw; j++) {
      addC[(i + (8 - j)) % 8] = -1.0 * std::numeric_limits<double>::max();
    }
  }
}

// _____________________________________________________________________________
void GridGraph::topoBlockPenalty(GridNode* n, CombNode* origNode, CombEdge* e,
                                 double* addC) {
  CombEdge* outgoing[8];
  getSettledOutgoingEdges(n, outgoing);

  // topological blocking
  for (size_t i = 0; i < 8; i++) {
    if (!outgoing[i]) continue;

    for (size_t j = i + 1; j < i + 8; j++) {
      if (!outgoing[j % 8]) continue;
      if (outgoing[j % 8] == outgoing[i]) break;

      int da = origNode->pl().getEdgeOrdering().dist(outgoing[i], e);
      int db = origNode->pl().getEdgeOrdering().dist(outgoing[j % 8], e);

      if (db < da) {
        // edge does not lie in this segment, block it!
        for (size_t x = i + 1; x < j; x++) {
          addC[x % 8] = -1.0 * std::numeric_limits<double>::max();
        }
      }
    }
  }
}

// _____________________________________________________________________________
void GridGraph::outDegDeviationPenalty(GridNode* n, CombNode* origNode,
                                       CombEdge* e, double* addC) {
  double degA = util::geo::angBetween(
      *origNode->pl().getParent()->pl().getGeom(),
      *e->getOtherNd(origNode)->pl().getParent()->pl().getGeom());

  int deg = -degA * (180.0 / M_PI);
  if (deg < 0) deg += 360;

  deg = (deg + 90) % 360;

  std::cerr << "deg is " << deg << "(degA: " << degA << ")" << std::endl;

  for (int i = 0; i < 8; i++) {
    if (addC[i] < -1) continue;
    double diff = std::min<int>(abs(deg - (45 * i)), 360 - abs(deg - (45 * i)));

    std::cerr << "diff @" << i << ": " << diff << std::endl;

    double multiplier = .1;
    addC[i] += multiplier * diff;
  }
}

// _____________________________________________________________________________
void GridGraph::addCostVector(GridNode* n, double addC[8], double* invCost) {
  std::cerr << "Adding cost vector ";
  for (size_t i = 0; i < 8; i++) {
    std::cerr << addC[i] << ",";
  }
  std::cerr << std::endl;

  auto xy = getNodeCoords(n);
  size_t x = xy.first;
  size_t y = xy.second;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    auto neigh = getNeighbor(x, y, i);

    if (!neigh) continue;

    auto oPort = neigh->pl().getPort((i + 4) % 8);

    if (addC[i] < -1) {
      if (getEdg(port, oPort)->pl().closed()) {
        // already closed, so dont remove closedness in inv costs
        invCost[i] = 0;
      } else {
        // close edges
        getEdg(port, oPort)->pl().close();
        getEdg(oPort, port)->pl().close();

        // close the other node to avoid "stealing" the edge
        // IMPORTANT: because we check if this edge is already closed
        // above, it is impossible for this node to be already closed -
        // then the edge would also be close, too. So we can close this node
        // here without danger of re-opening an already closed node later on.
        closeNode(getNeighbor(x, y, i));

        invCost[i] = addC[i];
      }
    } else {
      getEdg(port, oPort)
          ->pl()
          .setCost(getEdg(port, oPort)->pl().rawCost() + addC[i]);
      getEdg(oPort, port)
          ->pl()
          .setCost(getEdg(oPort, port)->pl().rawCost() + addC[i]);
      invCost[i] = addC[i];
    }
  }
}

// _____________________________________________________________________________
void GridGraph::removeCostVector(GridNode* n, double addC[8]) {
  std::cerr << "Removing cost vector ";
  for (size_t i = 0; i < 8; i++) {
    std::cerr << addC[i] << ",";
  }
  std::cerr << std::endl;
  auto xy = getNodeCoords(n);
  size_t x = xy.first;
  size_t y = xy.second;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    auto neigh = getNeighbor(x, y, i);

    if (!neigh) continue;

    auto oPort = neigh->pl().getPort((i + 4) % 8);

    if (addC[i] < -1) {
      getEdg(port, oPort)->pl().open();
      getEdg(oPort, port)->pl().open();
      openNode(getNeighbor(x, y, i));
    } else {
      getEdg(port, oPort)
          ->pl()
          .setCost(

              getEdg(port, oPort)->pl().rawCost() - addC[i]);
      getEdg(oPort, port)
          ->pl()
          .setCost(

              getEdg(oPort, port)->pl().rawCost() - addC[i]);
    }
  }
}

// _____________________________________________________________________________
std::set<util::graph::Edge<octi::combgraph::CombNodePL,
                           octi::combgraph::CombEdgePL>*>
GridGraph::getResEdges(GridNode* n) const {
  std::set<util::graph::Edge<octi::combgraph::CombNodePL,
                             octi::combgraph::CombEdgePL>*>
      ret;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    for (auto e : port->getAdjList()) {
      ret.insert(e->pl().getResEdges().begin(), e->pl().getResEdges().end());
    }
  }

  return ret;
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
        auto e = getEdg(port, neigh->pl().getPort((i + 4) % 8));
        auto f = getEdg(neigh->pl().getPort((i + 4) % 8), port);

        if (i % 4 == 0) {
          e->pl().setCost(_c.verticalPen);
          f->pl().setCost(_c.verticalPen);
        }
        if ((i + 2) % 4 == 0) {
          e->pl().setCost(_c.horizontalPen);
          f->pl().setCost(_c.horizontalPen);
        }
        if (i % 2) {
          e->pl().setCost(_c.diagonalPen);
          f->pl().setCost(_c.diagonalPen);
        }
      }
    }
  }
}

// _____________________________________________________________________________
std::priority_queue<Candidate> GridGraph::getNearestCandidatesFor(
    const util::geo::DPoint& p, double maxD) const {
  std::priority_queue<Candidate> ret;
  std::set<GridNode*> neigh;
  util::geo::DBox b(util::geo::DPoint(p.getX() - maxD, p.getY() - maxD),
                    util::geo::DPoint(p.getX() + maxD, p.getY() + maxD));
  _grid.get(b, &neigh);

  for (auto n : neigh) {
    if (n->pl().isClosed()) continue;
    double d = util::geo::dist(*n->pl().getGeom(), p);
    if (d < maxD) {
      ret.push(Candidate(n, d));
    }
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
  if (xa == xb && ya == yb) return 0;
  size_t minHops = std::max(labs(xb - xa), labs(yb - ya));

  size_t edgeCost =
      minHops *
      (std::min(_c.verticalPen, std::min(_c.horizontalPen, _c.diagonalPen)));
  size_t hopCost = (minHops - 1) * (_c.p_45 - _c.p_135);

  return edgeCost + hopCost;
}

// _____________________________________________________________________________
void GridGraph::openNode(GridNode* n) {
  if (!n->pl().isClosed()) return;
  auto xy = getNodeCoords(n);
  size_t x = xy.first;
  size_t y = xy.second;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    auto neigh = getNeighbor(x, y, i);
    if (!neigh || !port || neigh->pl().isClosed()) continue;
    auto e = getEdg(port, neigh->pl().getPort((i + 4) % 8));
    auto f = getEdg(neigh->pl().getPort((i + 4) % 8), port);
    if (e->pl().getResEdges().size() == 0) {
      e->pl().open();
      f->pl().open();
    }
  }

  n->pl().setClosed(false);
}

// _____________________________________________________________________________
void GridGraph::closeNode(GridNode* n) {
  if (n->pl().isClosed()) return;
  auto xy = getNodeCoords(n);
  size_t x = xy.first;
  size_t y = xy.second;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    auto neigh = getNeighbor(x, y, i);
    if (!neigh || !port) continue;
    auto e = getEdg(port, neigh->pl().getPort((i + 4) % 8));
    auto f = getEdg(neigh->pl().getPort((i + 4) % 8), port);

    e->pl().close();
    f->pl().close();
  }

  n->pl().setClosed(true);
}

// _____________________________________________________________________________
void GridGraph::openNodeSink(GridNode* n, double cost) {
  for (size_t i = 0; i < 8; i++) {
    getEdg(n->pl().getPort(i), n)->pl().setCost(cost);
    getEdg(n, n->pl().getPort(i))->pl().setCost(cost);
  }
}

// _____________________________________________________________________________
void GridGraph::closeNodeSink(GridNode* n) {
  for (size_t i = 0; i < 8; i++) {
    getEdg(n->pl().getPort(i), n)->pl().setCost(INF);
    getEdg(n, n->pl().getPort(i))->pl().setCost(INF);
  }
}

// _____________________________________________________________________________
GridNode* GridGraph::getGridNodeFrom(CombNode* n, double maxDis) {
  if (!isSettled(n)) {
    auto cands = getNearestCandidatesFor(*n->pl().getGeom(), maxDis);
    std::cerr << cands.size() << std::endl;

    while (!cands.empty()) {
      if (!cands.top().n->pl().isClosed()) {
        return cands.top().n;
      }
      cands.pop();
    }
    return 0;
  }
  return _settled[n];
}

// _____________________________________________________________________________
std::set<GridNode*> GridGraph::getGridNodesTo(CombNode* n, double maxDis) {
  std::set<GridNode*> tos;
  if (!isSettled(n)) {
    auto cands = getNearestCandidatesFor(*n->pl().getGeom(), maxDis);

    while (!cands.empty()) {
      if (!cands.top().n->pl().isClosed()) {
        tos.insert(cands.top().n);
      }
      cands.pop();
    }
  } else {
    tos.insert(_settled.find(n)->second);
  }

  return tos;
}

// _____________________________________________________________________________
void GridGraph::settleGridNode(GridNode* n, CombNode* cn) { _settled[cn] = n; }

// _____________________________________________________________________________
bool GridGraph::isSettled(CombNode* cn) {
  return _settled.find(cn) != _settled.end();
}
