// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/GridGraph.h"
#include "util/graph/Node.h"

using namespace octi::gridgraph;

// _____________________________________________________________________________
GridGraph::GridGraph(const util::geo::Box& bbox, double cellSize)
    : _bbox(bbox), _grid(cellSize, cellSize, bbox), Graph<NodePL, EdgePL>(true) {

  _c = {
    0, // cost for 0 degree node traversal
    90, // cost for 45 degree node traversal
    30, // cost for 90 degree node traversal
    10, // cost for 135 degree node traversal
    3, // cost for vertical edge
    3, // cost for horizontal edge
    3 // cost for diagonal edge
  };

  assert(_c.p_0 < _c.p_135);
  assert(_c.p_135 < _c.p_90);
  assert(_c.p_90 < _c.p_45);

  double c_0 = _c.p_45 - _c.p_135;
  double c_135 = _c.p_45;
  double c_90 = _c.p_45 - _c.p_135 + _c.p_90;

  double directNodeCost = 88888;

  // write nodes
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      double xPos = bbox.min_corner().get<0>() + x * cellSize;
      double yPos = bbox.min_corner().get<1>() + y * cellSize;
      Node<NodePL, EdgePL>* n = new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos)));
      _grid.add(x, y, n);
      addNode(n);
      n->pl().setParent(n);

      Node<NodePL, EdgePL>* n1 = new Node<NodePL, EdgePL>(NodePL(Point(xPos - cellSize/10.0, yPos)));
      n1->pl().setParent(n);
      addNode(n1);
      n->pl().setWestPort(n1);
      addEdge(n, n1, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n1->pl().getGeom()), directNodeCost));
      addEdge(n1, n, EdgePL(util::geo::PolyLine(*n1->pl().getGeom(), *n->pl().getGeom()), directNodeCost));

      Node<NodePL, EdgePL>* n2 = new Node<NodePL, EdgePL>(NodePL(Point(xPos - cellSize/10.0, yPos - cellSize/10.0)));
      addNode(n2);
      n2->pl().setParent(n);
      n->pl().setSouthWestPort(n2);
      addEdge(n, n2, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n2->pl().getGeom()), directNodeCost));
      addEdge(n2, n, EdgePL(util::geo::PolyLine(*n2->pl().getGeom(), *n->pl().getGeom()), directNodeCost));

      Node<NodePL, EdgePL>* n3 = new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos - cellSize/10.0)));
      addNode(n3);
      n3->pl().setParent(n);
      n->pl().setSouthPort(n3);
      addEdge(n, n3, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n3->pl().getGeom()), directNodeCost));
      addEdge(n3, n, EdgePL(util::geo::PolyLine(*n3->pl().getGeom(), *n->pl().getGeom()), directNodeCost));

      Node<NodePL, EdgePL>* n4 = new Node<NodePL, EdgePL>(NodePL(Point(xPos + cellSize/10.0, yPos)));
      addNode(n4);
      n4->pl().setParent(n);
      n->pl().setEastPort(n4);
      addEdge(n, n4, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n4->pl().getGeom()), directNodeCost));
      addEdge(n4, n, EdgePL(util::geo::PolyLine(*n4->pl().getGeom(), *n->pl().getGeom()), directNodeCost));

      Node<NodePL, EdgePL>* n5 = new Node<NodePL, EdgePL>(NodePL(Point(xPos + cellSize/10.0, yPos + cellSize/10.0)));
      n5->pl().setParent(n);
      addNode(n5);
      n->pl().setNorthEastPort(n5);
      addEdge(n, n5, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n5->pl().getGeom()), directNodeCost));
      addEdge(n5, n, EdgePL(util::geo::PolyLine(*n5->pl().getGeom(), *n->pl().getGeom()), directNodeCost));

      Node<NodePL, EdgePL>* n6 = new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos + cellSize/10.0)));
      n6->pl().setParent(n);
      addNode(n6);
      n->pl().setNorthPort(n6);
      addEdge(n, n6, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n6->pl().getGeom()), directNodeCost));
      addEdge(n6, n, EdgePL(util::geo::PolyLine(*n6->pl().getGeom(), *n->pl().getGeom()), directNodeCost));

      Node<NodePL, EdgePL>* n7 = new Node<NodePL, EdgePL>(NodePL(Point(xPos - cellSize/10.0, yPos + cellSize/10.0)));
      n7->pl().setParent(n);
      addNode(n7);
      n->pl().setNorthWestPort(n7);
      addEdge(n, n7, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n7->pl().getGeom()), directNodeCost));
      addEdge(n7, n, EdgePL(util::geo::PolyLine(*n7->pl().getGeom(), *n->pl().getGeom()), directNodeCost));

      Node<NodePL, EdgePL>* n8 = new Node<NodePL, EdgePL>(NodePL(Point(xPos + cellSize/10.0, yPos - cellSize/10.0)));
      n8->pl().setParent(n);
      addNode(n8);
      n->pl().setSouthEastPort(n8);
      addEdge(n, n8, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n8->pl().getGeom()), directNodeCost));
      addEdge(n8, n, EdgePL(util::geo::PolyLine(*n8->pl().getGeom(), *n->pl().getGeom()), directNodeCost));


      // in-node connections
      for (size_t i = 0; i < 8; i++) {
        for (size_t j = i + 1; j < 8; j++) {
          int d = (int)(i) - (int)(j);
          size_t deg = abs((((d + 4) % 8) + 8) % 8 - 4);
          double pen = c_0;

          if (deg == 1) continue;
          if (deg == 2) pen = c_90;
          if (deg == 3) pen = c_135;
          addEdge(n->pl().getPort(i), n->pl().getPort(j), EdgePL(util::geo::PolyLine(*n->pl().getPort(i)->pl().getGeom(), *n->pl().getPort(j)->pl().getGeom()), pen));
          addEdge(n->pl().getPort(j), n->pl().getPort(i), EdgePL(util::geo::PolyLine(*n->pl().getPort(j)->pl().getGeom(), *n->pl().getPort(i)->pl().getGeom()), pen));
        }
      }

    }
  }

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      Node<NodePL, EdgePL>* center = getNode(x, y);
      if (center == 0) continue;

      Node<NodePL, EdgePL>* from = center->pl().getNorthPort();
      Node<NodePL, EdgePL>* toN = getNeighbor(x, y, 0);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getSouthPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }

      from = center->pl().getSouthPort();
      toN = getNeighbor(x, y, 4);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getNorthPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }

      from = center->pl().getWestPort();
      toN = getNeighbor(x, y, 6);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getEastPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }

      from = center->pl().getEastPort();
      toN = getNeighbor(x, y, 2);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getWestPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }

      from = center->pl().getSouthWestPort();
      toN = getNeighbor(x, y, 5);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getNorthEastPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }

      from = center->pl().getSouthEastPort();
      toN = getNeighbor(x, y, 3);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getNorthWestPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }

      from = center->pl().getNorthWestPort();
      toN = getNeighbor(x, y, 7);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getSouthEastPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }

      from = center->pl().getNorthEastPort();
      toN = getNeighbor(x, y, 1);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getSouthWestPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 0));
      }
    }
  }
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
Node<NodePL, EdgePL>* GridGraph::getNeighbor(size_t cx, size_t cy, size_t i) const {
  int8_t x = 1;
  if (i % 4 == 0) x = 0;
  if (i > 4) x = -1;

  int8_t y = 1;
  if (i == 2 || i == 6) y = 0;
  if (i == 3 || i == 4 || i == 5) y = -1;

  return getNode(cx + x, cy + y);
}

// _____________________________________________________________________________
void GridGraph::balance() {
  writeInitialCosts();

  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      auto n = getNode(x, y);

      size_t numRoutes = getRoutes(n).size();

      for (size_t i = 0; i < 8; i++) {
        auto port = n->pl().getPort(i);
        auto neigh = getNeighbor(x, y, i);
        if (!neigh || !port) continue;
        auto e = getEdge(port, neigh->pl().getPort((i + 4) % 8));
        e->pl().setCost(e->pl().cost() + 99999 * numRoutes);

        if (i == 1 || i == 3) {
          size_t numRoutesEdge = e->pl().getRoutes().size();
          auto a = getNeighbor(x, y, (i + 7) % 8);
          auto b = getNeighbor(x, y, (i + 1) % 8);

          if (a && b) {
            size_t port = (i + 2) % 8;
            auto e = getEdge(a->pl().getPort(port), b->pl().getPort((port + 4) % 8));
            auto f = getEdge(b->pl().getPort((port + 4) % 8), a->pl().getPort(port));

            e->pl().setCost(e->pl().cost() + 99999 * numRoutesEdge);
            f->pl().setCost(f->pl().cost() + 99999 * numRoutesEdge);
          }
        }
      }

      // if some outgoing edge is taken, dont put new edge next to it
      for (size_t i = 0; i < 8; i++) {
        auto port = n->pl().getPort(i);
        auto neigh = getNeighbor(x, y, i);
        if (!neigh || !port) continue;
        auto e = getEdge(port, neigh->pl().getPort((i + 4) % 8));

        if (e->pl().getRoutes().size()) {
          if (getNeighbor(x, y, (i+1) % 8)) {
            auto e1 = getEdge(n->pl().getPort((i+1) % 8), getNeighbor(x, y, (i+1) % 8)->pl().getPort((i + 1 + 4) % 8));
            e1->pl().setCost(e1->pl().cost() + 150);
          }

          if (getNeighbor(x, y, (i+7) % 8)) {
            auto e2 = getEdge(n->pl().getPort((i+7) % 8), getNeighbor(x, y, (i+7) % 8)->pl().getPort((i + 7 + 4) % 8));
            e2->pl().setCost(e2->pl().cost() + 150);
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
std::set<std::string> GridGraph::getRoutes(Node<NodePL, EdgePL>* n) const {
  std::set<std::string> ret;

  for (size_t i = 0; i < 8; i++) {
    auto port = n->pl().getPort(i);
    for (auto e : port->getAdjList()) {
      ret.insert(e->pl().getRoutes().begin(), e->pl().getRoutes().end());
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

        if (i % 4 == 0) e->pl().setCost(_c.verticalCost);
        if ((i + 2) % 4 == 0) e->pl().setCost(_c.horizontalCost);
        if (i % 2) e->pl().setCost(_c.diagonalCost);
      }
    }
  }
}

// _____________________________________________________________________________
std::priority_queue<Candidate> GridGraph::getNearestCandidatesFor(const util::geo::Point& p, double maxD) const {
  std::priority_queue<Candidate> ret;
  std::set<Node<NodePL, EdgePL>*> neigh;
  util::geo::Box b(util::geo::Point(p.get<0>() - maxD, p.get<1>() - maxD), util::geo::Point(p.get<0>() + maxD, p.get<1>() + maxD));
  _grid.get(b, &neigh);

  for (auto n : neigh) {
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
