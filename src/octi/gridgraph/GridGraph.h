// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/Graph.h"
#include "octi/gridgraph/NodePL.h"
#include "octi/gridgraph/EdgePL.h"

using util::graph::Graph;
using util::graph::Node;
using util::geo::Grid;
using util::geo::Point;

namespace octi {
namespace gridgraph {

class GridGraph : public Graph<NodePL, EdgePL> {
 public:
  GridGraph(const util::geo::Box& bbox, double cellSize);

  Node<NodePL, EdgePL>* getNode(size_t x, size_t y) const;

 private:
  util::geo::Box _bbox;

  Grid<Node<NodePL, EdgePL>*, Point> _grid;

  Node<NodePL, EdgePL>* getN(Node<NodePL, EdgePL>* n) const;
  Node<NodePL, EdgePL>* getNE(Node<NodePL, EdgePL>* n) const;
  Node<NodePL, EdgePL>* getE(Node<NodePL, EdgePL>* n) const;
  Node<NodePL, EdgePL>* getSE(Node<NodePL, EdgePL>* n) const;
  Node<NodePL, EdgePL>* getS(Node<NodePL, EdgePL>* n) const;
  Node<NodePL, EdgePL>* getSW(Node<NodePL, EdgePL>* n) const;
  Node<NodePL, EdgePL>* getW(Node<NodePL, EdgePL>* n) const;
  Node<NodePL, EdgePL>* getNW(Node<NodePL, EdgePL>* n) const;

};

}
}
