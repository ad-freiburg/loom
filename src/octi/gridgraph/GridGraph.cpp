// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include "octi/gridgraph/GridGraph.h"
#include "util/graph/Node.h"

using namespace octi::gridgraph;

double INF = std::numeric_limits<double>::infinity();

// _____________________________________________________________________________
GridGraph::GridGraph(const util::geo::Box& bbox, double cellSize,
                     const Penalties& pens)
    : _bbox(bbox),
      _c(pens),
      _grid(cellSize, cellSize, bbox),
      Graph<NodePL, EdgePL>(true) {

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
      double xPos = bbox.min_corner().get<0>() + x * cellSize;
      double yPos = bbox.min_corner().get<1>() + y * cellSize;
      Node<NodePL, EdgePL>* n =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos)));
      _grid.add(x, y, n);
      n->pl().setXY(x, y);
      addNode(n);
      n->pl().setParent(n);

      Node<NodePL, EdgePL>* n1 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos - spacer, yPos)));
      n1->pl().setParent(n);
      addNode(n1);
      n->pl().setPort(6, n1);
      addEdge(n, n1, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n1->pl().getGeom()),
                            INF, true, false));
      addEdge(n1, n, EdgePL(util::geo::PolyLine(*n1->pl().getGeom(),
                                                *n->pl().getGeom()),
                            INF, true, false));

      Node<NodePL, EdgePL>* n2 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos - spacer, yPos - spacer)));
      addNode(n2);
      n2->pl().setParent(n);
      n->pl().setPort(5, n2);
      addEdge(n, n2, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n2->pl().getGeom()),
                            INF, true, false));
      addEdge(n2, n, EdgePL(util::geo::PolyLine(*n2->pl().getGeom(),
                                                *n->pl().getGeom()),
                            INF, true, false));

      Node<NodePL, EdgePL>* n3 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos - spacer)));
      addNode(n3);
      n3->pl().setParent(n);
      n->pl().setPort(4, n3);
      addEdge(n, n3, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n3->pl().getGeom()),
                            INF, true, false));
      addEdge(n3, n, EdgePL(util::geo::PolyLine(*n3->pl().getGeom(),
                                                *n->pl().getGeom()),
                            INF, true, false));

      Node<NodePL, EdgePL>* n4 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos + spacer, yPos)));
      addNode(n4);
      n4->pl().setParent(n);
      n->pl().setPort(2, n4);
      addEdge(n, n4, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n4->pl().getGeom()),
                            INF, true, false));
      addEdge(n4, n, EdgePL(util::geo::PolyLine(*n4->pl().getGeom(),
                                                *n->pl().getGeom()),
                            INF, true, false));

      Node<NodePL, EdgePL>* n5 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos + spacer, yPos + spacer)));
      n5->pl().setParent(n);
      addNode(n5);
      n->pl().setPort(1, n5);
      addEdge(n, n5, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n5->pl().getGeom()),
                            INF, true, false));
      addEdge(n5, n, EdgePL(util::geo::PolyLine(*n5->pl().getGeom(),
                                                *n->pl().getGeom()),
                            INF, true, false));

      Node<NodePL, EdgePL>* n6 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos + spacer)));
      n6->pl().setParent(n);
      addNode(n6);
      n->pl().setPort(0, n6);
      addEdge(n, n6, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n6->pl().getGeom()),
                            INF, true, false));
      addEdge(n6, n, EdgePL(util::geo::PolyLine(*n6->pl().getGeom(),
                                                *n->pl().getGeom()),
                            INF, true, false));

      Node<NodePL, EdgePL>* n7 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos - spacer, yPos + spacer)));
      n7->pl().setParent(n);
      addNode(n7);
      n->pl().setPort(7, n7);
      addEdge(n, n7, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n7->pl().getGeom()),
                            INF, true, false));
      addEdge(n7, n, EdgePL(util::geo::PolyLine(*n7->pl().getGeom(),
                                                *n->pl().getGeom()),
                            INF, true, false));
      Node<NodePL, EdgePL>* n8 =
          new Node<NodePL, EdgePL>(NodePL(Point(xPos + spacer, yPos - spacer)));
      n8->pl().setParent(n);
      addNode(n8);
      n->pl().setPort(3, n8);
      addEdge(n, n8, EdgePL(util::geo::PolyLine(*n->pl().getGeom(),
                                                *n8->pl().getGeom()),
                            INF, true, false));
      addEdge(n8, n, EdgePL(util::geo::PolyLine(*n8->pl().getGeom(),
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
          addEdge(
              n->pl().getPort(i), n->pl().getPort(j),
              EdgePL(util::geo::PolyLine(*n->pl().getPort(i)->pl().getGeom(),
                                         *n->pl().getPort(j)->pl().getGeom()),
                     pen, true));
          addEdge(
              n->pl().getPort(j), n->pl().getPort(i),
              EdgePL(util::geo::PolyLine(*n->pl().getPort(j)->pl().getGeom(),
                                         *n->pl().getPort(i)->pl().getGeom()),
                     pen, true));
        }
      }
    }
  }

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      Node<NodePL, EdgePL>* center = getNode(x, y);
      if (center == 0) continue;

      for (size_t p = 0; p < 8; p++) {
        Node<NodePL, EdgePL>* from = center->pl().getPort(p);
        Node<NodePL, EdgePL>* toN = getNeighbor(x, y, p);
        if (from != 0 && toN != 0) {
          Node<NodePL, EdgePL>* to = toN->pl().getPort((p + 4) % 8);
          addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(),
                                                       *to->pl().getGeom()),
                                   0, false));
        }
      }
    }
  }

  writeInitialCosts();
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getNode(size_t x, size_t y) const {
  if (x >= _grid.getXWidth() || y >= _grid.getYHeight()) return 0;
  std::set<Node<NodePL, EdgePL>*> r;
  _grid.get(x, y, &r);

  if (r.size()) return *r.begin();
  return 0;
}

// _____________________________________________________________________________
std::pair<size_t, size_t> GridGraph::getNodeCoords(GridNode* n) const {
  return std::pair<size_t, size_t>(n->pl().getX(), n->pl().getY());
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getNeighbor(size_t cx, size_t cy,
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
    if (getEdge(a->pl().getPort(dir), b->pl().getPort((dir + 4) % 8))) {
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

  double nearPenalty = 10;

  // std::cerr << x << ", " << y << std::endl;

  if (dir == 0) {
    auto a = getNEdge(getNode(x - 1, y), getNode(x - 1, y + 1));
    auto b = getNEdge(getNode(x + 1, y), getNode(x + 1, y + 1));

    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty);
    }
  } else if (dir == 4) {
    auto a = getNEdge(getNode(x - 1, y), getNode(x - 1, y - 1));
    auto b = getNEdge(getNode(x + 1, y), getNode(x + 1, y - 1));

    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty);
    }
  }
  if (dir == 2) {
    auto a = getNEdge(getNode(x, y + 1), getNode(x + 1, y + 1));
    auto b = getNEdge(getNode(x, y - 1), getNode(x + 1, y - 1));

    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty);
    }
  } else if (dir == 6) {
    auto a = getNEdge(getNode(x, y + 1), getNode(x - 1, y + 1));
    auto b = getNEdge(getNode(x, y - 1), getNode(x - 1, y - 1));

    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty);
    }

    // diagonal
  } else if (dir == 1) {
    auto a = getNEdge(getNode(x - 1, y), getNode(x, y + 1));
    auto b = getNEdge(getNode(x, y + 1), getNode(x + 1, y + 2));
    auto c = getNEdge(getNode(x, y - 1), getNode(x + 1, y));
    auto d = getNEdge(getNode(x + 1, y), getNode(x + 2, y + 1));
    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (c) {
      c->pl().setCost(c->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(c)->pl().setCost(getOtherEdge(c)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (d) {
      d->pl().setCost(d->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(d)->pl().setCost(getOtherEdge(d)->pl().rawCost() +
                                    nearPenalty / 2);
    }
  } else if (dir == 5) {
    auto a = getNEdge(getNode(x - 1, y), getNode(x, y + 1));
    auto b = getNEdge(getNode(x - 2, y - 1), getNode(x - 1, y));
    auto c = getNEdge(getNode(x, y - 1), getNode(x + 1, y));
    auto d = getNEdge(getNode(x - 1, y - 2), getNode(x, y - 1));
    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (c) {
      c->pl().setCost(c->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(c)->pl().setCost(getOtherEdge(c)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (d) {
      d->pl().setCost(d->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(d)->pl().setCost(getOtherEdge(d)->pl().rawCost() +
                                    nearPenalty / 2);
    }
  } else if (dir == 3) {
    auto a = getNEdge(getNode(x - 1, y), getNode(x, y - 1));
    auto b = getNEdge(getNode(x, y - 1), getNode(x + 1, y - 2));
    auto c = getNEdge(getNode(x, y + 1), getNode(x + 1, y));
    auto d = getNEdge(getNode(x + 1, y), getNode(x + 2, y - 1));
    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (c) {
      c->pl().setCost(c->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(c)->pl().setCost(getOtherEdge(c)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (d) {
      d->pl().setCost(d->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(d)->pl().setCost(getOtherEdge(d)->pl().rawCost() +
                                    nearPenalty / 2);
    }
  } else if (dir == 7) {
    auto a = getNEdge(getNode(x, y - 1), getNode(x - 1, y));
    auto b = getNEdge(getNode(x - 1, y), getNode(x - 2, y + 1));
    auto c = getNEdge(getNode(x + 1, y), getNode(x, y + 1));
    auto d = getNEdge(getNode(x, y + 1), getNode(x - 1, y + 2));
    if (a) {
      a->pl().setCost(a->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(a)->pl().setCost(getOtherEdge(a)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (b) {
      b->pl().setCost(b->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(b)->pl().setCost(getOtherEdge(b)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (c) {
      c->pl().setCost(c->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(c)->pl().setCost(getOtherEdge(c)->pl().rawCost() +
                                    nearPenalty / 2);
    }
    if (d) {
      d->pl().setCost(d->pl().rawCost() + nearPenalty / 2);
      getOtherEdge(d)->pl().setCost(getOtherEdge(d)->pl().rawCost() +
                                    nearPenalty / 2);
    }
  }

  // std::cerr << dir << std::endl;
}

// _____________________________________________________________________________
GridEdge* GridGraph::getOtherEdge(GridEdge* e) {
  return getEdge(e->getTo(), e->getFrom());
}

// _____________________________________________________________________________
GridEdge* GridGraph::getNEdge(GridNode* a, GridNode* b) {
  for (size_t i = 0; i < 8; i++) {
    if (a && b && a->pl().getPort(i) && b->pl().getPort((i + 4) % 8)) {
      auto e = getEdge(a->pl().getPort(i), b->pl().getPort((i + 4) % 8));
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
        getEdge(port, neigh->pl().getPort((i + 4) % 8))
                ->pl()
                .getResEdges()
                .size() > 0) {
      outgoing[i] = *getEdge(port, neigh->pl().getPort((i + 4) % 8))
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

  if (origNode->pl().getParent()->pl().getStops().size())
    std::cerr << std::endl
              << std::endl
              << "Checking station "
              << origNode->pl().getParent()->pl().getStops().front().name
              << std::endl;
  else
    std::cerr << std::endl
              << std::endl
              << "Checking node " << origNode << std::endl;

  int origEdgeNumber = origNode->getAdjList().size();
  size_t optimDistance = (8 / origEdgeNumber) - 1;

  if (!origNode->pl().hasOrderedEdge(e)) {
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

  std::cerr << "Ordered edges in orig node: ";
  for (auto e : origNode->pl().getOrderedEdges()) {
    if (e.first->getOtherNode(origNode)
            ->pl()
            .getParent()
            ->pl()
            .getStops()
            .size()) {
      std::cerr << e.first << "(to "
                << e.first->getOtherNode(origNode)
                       ->pl()
                       .getParent()
                       ->pl()
                       .getStops()
                       .front()
                       .name
                << "), ";
    } else {
      std::cerr << e.first << ", ";
    }
  }
  std::cerr << std::endl;

  for (size_t i = 0; i < 8; i++) {
    if (!outgoing[i]) continue;

    // this is the number of edges that will occur between the currently checked
    // edge and the inserted edge,
    // positive if in clockwise direction, negative if in counter-clockwise
    // direction
    int32_t d = origNode->pl().distBetween(outgoing[i], e) - 1;

    // dd and ddd are the optimal distances between outgoing[i] and e, based on
    // exit(0);
    // the total number
    // of edges in this node
    int dd = ((((d + 1) + d) % 8) * optimDistance) % 8;
    int ddd = (6 - dd) % 8;

    std::cerr << "Distance between the inserted edge (" << e << ") and edge at "
              << i << " (" << outgoing[i] << ") is " << d
              << ", optim distance between them is +" << dd << " and -" << ddd
              << std::endl;

    double pen = _c.p_45 * 2 - 1;

    for (int j = 1; j <= dd + 1; j++) {
      if (addC[(i + j) % 8] < -1) continue;
      addC[(i + j) % 8] += pen * (1.0 - (j - 1.0) / (dd));
    }

    for (int j = 1; j <= ddd + 1; j++) {
      if (addC[(i + (8 - j)) % 8] < -1) continue;
      addC[(i + (8 - j)) % 8] += pen * (1.0 - (j - 1.0) / (ddd));
    }

    // negative cost here means that the edge is going to be closed
    addC[i] = -1.0 * std::numeric_limits<double>::max();

    for (int j = 1; j <= d; j++) {
      addC[(i + j) % 8] = -1.0 * std::numeric_limits<double>::max();
    }

    for (int j = 1; j <= (origEdgeNumber / 2) - d; j++) {
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

      int da = origNode->pl().distBetween(outgoing[i], e);
      int db = origNode->pl().distBetween(outgoing[j % 8], e);

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
void GridGraph::outDegDeviationPenalty(GridNode* n, CombNode* origNode, CombEdge* e,
                            double* addC) {
  double degA = util::geo::angBetween(
      *origNode->pl().getParent()->pl().getGeom(),
      *e->getOtherNode(origNode)->pl().getParent()->pl().getGeom());

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
  for (size_t i = 0 ; i < 8; i++) {
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
      if (getEdge(port, oPort)->pl().closed()) {
        // already closed, so dont remove closedness in inv costs
        invCost[i] = 0;
      } else {
        getEdge(port, oPort)->pl().close();
        getEdge(oPort, port)->pl().close();
        invCost[i] = addC[i];
      }
    } else {
      getEdge(port, oPort)
          ->pl()
          .setCost(
              getEdge(port, oPort)->pl().rawCost() + addC[i]);
      getEdge(oPort, port)
          ->pl()
          .setCost(
              getEdge(oPort, port)->pl().rawCost() + addC[i]);
      invCost[i] = addC[i];
    }
  }
}

// _____________________________________________________________________________
void GridGraph::removeCostVector(GridNode* n, double addC[8]) {
  auto xy = getNodeCoords(n);
  size_t x = xy.first;
  size_t y = xy.second;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    auto neigh = getNeighbor(x, y, i);

    if (!neigh) continue;

    auto oPort = neigh->pl().getPort((i + 4) % 8);

    if (addC[i] < -1) {
      getEdge(port, oPort)->pl().open();
      getEdge(oPort, port)->pl().open();
    } else {
      getEdge(port, oPort)
          ->pl()
          .setCost(

              getEdge(port, oPort)->pl().rawCost() - addC[i]);
      getEdge(oPort, port)
          ->pl()
          .setCost(

              getEdge(oPort, port)->pl().rawCost() - addC[i]);
    }
  }
}

// _____________________________________________________________________________
std::set<util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>*>
GridGraph::getResEdges(Node<NodePL, EdgePL>* n) const {
  std::set<util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>*>
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
        auto e = getEdge(port, neigh->pl().getPort((i + 4) % 8));
        auto f = getEdge(neigh->pl().getPort((i + 4) % 8), port);

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
    const util::geo::Point& p, double maxD) const {
  std::priority_queue<Candidate> ret;
  std::set<Node<NodePL, EdgePL>*> neigh;
  util::geo::Box b(util::geo::Point(p.get<0>() - maxD, p.get<1>() - maxD),
                   util::geo::Point(p.get<0>() + maxD, p.get<1>() + maxD));
  _grid.get(b, &neigh);

  for (auto n : neigh) {
    if (n->pl().isClosed()) continue;
    double d = util::geo::dist(n->pl().getGeom(), p);
    if (d < maxD) {
      ret.push(Candidate(n, d));
    }
  }

  return ret;
}

// _____________________________________________________________________________
const Grid<Node<NodePL, EdgePL>*, Point>& GridGraph::getGrid() const {
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
    if (!neigh || !port) continue;
    auto e = getEdge(port, neigh->pl().getPort((i + 4) % 8));
    auto f = getEdge(neigh->pl().getPort((i + 4) % 8), port);
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
    auto e = getEdge(port, neigh->pl().getPort((i + 4) % 8));
    auto f = getEdge(neigh->pl().getPort((i + 4) % 8), port);

    e->pl().close();
    f->pl().close();
  }

  n->pl().setClosed(true);
}

// _____________________________________________________________________________
void GridGraph::openNodeSink(GridNode* n, double cost) {
  for (size_t i = 0; i < 8; i++) {
    getEdge(n->pl().getPort(i), n)->pl().setCost(cost);
    getEdge(n, n->pl().getPort(i))->pl().setCost(cost);
  }
}

// _____________________________________________________________________________
void GridGraph::closeNodeSink(GridNode* n) {
  for (size_t i = 0; i < 8; i++) {
    getEdge(n->pl().getPort(i), n)->pl().setCost(INF);
    getEdge(n, n->pl().getPort(i))->pl().setCost(INF);
  }
}
