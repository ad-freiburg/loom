// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/gridgraph/GridGraph.h"
#include "util/graph/Node.h"

using namespace octi::gridgraph;

// _____________________________________________________________________________
GridGraph::GridGraph(const util::geo::Box& bbox, double cellSize)
    : _bbox(bbox), _grid(cellSize, cellSize, bbox), Graph<NodePL, EdgePL>(true) {

  // write nodes
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      double xPos = x * cellSize;
      double yPos = y * cellSize;
      Node<NodePL, EdgePL>* n = new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos)));
      _grid.add(x, y, n);
      addNode(n);

      Node<NodePL, EdgePL>* n1 = new Node<NodePL, EdgePL>(NodePL(Point(xPos - cellSize/10.0, yPos)));
      addNode(n1);
      n->pl().setWestPort(n1);
      addEdge(n, n1, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n1->pl().getGeom()), 5000));
      addEdge(n1, n, EdgePL(util::geo::PolyLine(*n1->pl().getGeom(), *n->pl().getGeom()), 5000));

      Node<NodePL, EdgePL>* n2 = new Node<NodePL, EdgePL>(NodePL(Point(xPos - cellSize/10.0, yPos - cellSize/10.0)));
      addNode(n2);
      n->pl().setSouthWestPort(n2);
      addEdge(n, n2, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n2->pl().getGeom()), 5000));
      addEdge(n2, n, EdgePL(util::geo::PolyLine(*n2->pl().getGeom(), *n->pl().getGeom()), 5000));

      Node<NodePL, EdgePL>* n3 = new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos - cellSize/10.0)));
      addNode(n3);
      n->pl().setSouthPort(n3);
      addEdge(n, n3, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n3->pl().getGeom()), 5000));
      addEdge(n3, n, EdgePL(util::geo::PolyLine(*n3->pl().getGeom(), *n->pl().getGeom()), 5000));

      Node<NodePL, EdgePL>* n4 = new Node<NodePL, EdgePL>(NodePL(Point(xPos + cellSize/10.0, yPos)));
      addNode(n4);
      n->pl().setEastPort(n4);
      addEdge(n, n4, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n4->pl().getGeom()), 5000));
      addEdge(n4, n, EdgePL(util::geo::PolyLine(*n4->pl().getGeom(), *n->pl().getGeom()), 5000));

      Node<NodePL, EdgePL>* n5 = new Node<NodePL, EdgePL>(NodePL(Point(xPos + cellSize/10.0, yPos + cellSize/10.0)));
      addNode(n5);
      n->pl().setNorthEastPort(n5);
      addEdge(n, n5, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n5->pl().getGeom()), 5000));
      addEdge(n5, n, EdgePL(util::geo::PolyLine(*n5->pl().getGeom(), *n->pl().getGeom()), 5000));

      Node<NodePL, EdgePL>* n6 = new Node<NodePL, EdgePL>(NodePL(Point(xPos, yPos + cellSize/10.0)));
      addNode(n6);
      n->pl().setNorthPort(n6);
      addEdge(n, n6, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n6->pl().getGeom()), 5000));
      addEdge(n6, n, EdgePL(util::geo::PolyLine(*n6->pl().getGeom(), *n->pl().getGeom()), 5000));

      Node<NodePL, EdgePL>* n7 = new Node<NodePL, EdgePL>(NodePL(Point(xPos - cellSize/10.0, yPos + cellSize/10.0)));
      addNode(n7);
      n->pl().setNorthWestPort(n7);
      addEdge(n, n7, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n7->pl().getGeom()), 5000));
      addEdge(n7, n, EdgePL(util::geo::PolyLine(*n7->pl().getGeom(), *n->pl().getGeom()), 5000));

      Node<NodePL, EdgePL>* n8 = new Node<NodePL, EdgePL>(NodePL(Point(xPos + cellSize/10.0, yPos - cellSize/10.0)));
      addNode(n8);
      n->pl().setSouthEastPort(n8);
      addEdge(n, n8, EdgePL(util::geo::PolyLine(*n->pl().getGeom(), *n8->pl().getGeom()), 5000));
      addEdge(n8, n, EdgePL(util::geo::PolyLine(*n8->pl().getGeom(), *n->pl().getGeom()), 5000));


      // in-node connections
      for (size_t i = 0; i < 8; i++) {
        for (size_t j = i + 1; j < 8; j++) {
          int d = (int)(i) - (int)(j);
          double pen = 4 - fabs((((d + 4) % 8) + 8) % 8 - 4);
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
      Node<NodePL, EdgePL>* toN = getN(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getSouthPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
      }

      from = center->pl().getSouthPort();
      toN = getS(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getNorthPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
      }

      from = center->pl().getWestPort();
      toN = getW(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getEastPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
      }

      from = center->pl().getEastPort();
      toN = getE(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getWestPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
      }

      from = center->pl().getSouthWestPort();
      toN = getSW(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getNorthEastPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
      }

      from = center->pl().getSouthEastPort();
      toN = getSE(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getNorthWestPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
      }

      from = center->pl().getNorthWestPort();
      toN = getNW(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getSouthEastPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
      }

      from = center->pl().getNorthEastPort();
      toN = getNE(center);
      if (from != 0 && toN != 0) {
        Node<NodePL, EdgePL>* to = toN->pl().getSouthWestPort();
        addEdge(from, to, EdgePL(util::geo::PolyLine(*from->pl().getGeom(), *to->pl().getGeom()), 1));
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
Node<NodePL, EdgePL>* GridGraph::getNE(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x+1, y+1);
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getE(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x+1, y);
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getSE(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x+1, y-1);
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getS(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x, y-1);
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getSW(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x-1, y-1);
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getW(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x-1, y);
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getNW(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x-1, y+1);
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getN(Node<NodePL, EdgePL>* n) const {
  auto cells = _grid.getCells(n);
  size_t x = cells.begin()->first;
  size_t y = cells.begin()->second;

  return getNode(x, y+1);
}