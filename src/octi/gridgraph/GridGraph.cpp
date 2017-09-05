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
      Node<NodePL, EdgePL>* n = new Node<NodePL, EdgePL>(NodePL(Point(x, y)));
      _grid.add(x, y, n);
      addNode(n);
    }
  }

  std::set<Node<NodePL, EdgePL>*> neighbors;
  std::set<Node<NodePL, EdgePL>*> from;
  for (size_t x = 0; x < _grid.getXWidth(); x++) {
    for (size_t y = 0; y < _grid.getYHeight(); y++) {
      _grid.get(x, y, &from);
      Node<NodePL, EdgePL>* f = *from.begin();
      from.clear();
      _grid.getCellNeighbors(x, y, 1, 1, &neighbors);
      for (Node<NodePL, EdgePL>* n : neighbors) {
        if (n != f) {
          addEdge(f, n, EdgePL(util::geo::PolyLine(*f->pl().getGeom(), *n->pl().getGeom())));
        }
      }

      neighbors.clear();
    }
  }
}

// _____________________________________________________________________________
Node<NodePL, EdgePL>* GridGraph::getNode(size_t x, size_t y) const {
  std::set<Node<NodePL, EdgePL>*> r;
  _grid.get(x, y, &r);

  if (r.size()) return *r.begin();
  return 0;
}
