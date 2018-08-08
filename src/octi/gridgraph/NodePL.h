// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_NODEPL_H_
#define OCTI_GRIDGRAPH_NODEPL_H_

#include "octi/gridgraph/EdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace gridgraph {

class NodePL : util::geograph::GeoNodePL<double> {
 public:
  NodePL() {};
  NodePL(Point<double> pos);

  const Point<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  Node<NodePL, EdgePL>* getParent() const;
  void setParent(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getPort(size_t i) const;
  void setPort(size_t p, Node<NodePL, EdgePL>* n);

  void setXY(size_t x, size_t y);
  size_t getX() const;
  size_t getY() const;

  bool isClosed() const;
  void setClosed(bool c);


 private:
  Point<double> _pos;

  Node<NodePL, EdgePL>* _parent;

  Node<NodePL, EdgePL>* _ports[8];

  size_t _x, _y;
  bool _closed;
};
}
}

#endif  // OCTI_GRIDGRAPH_NODEPL_H_
