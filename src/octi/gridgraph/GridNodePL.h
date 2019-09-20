// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_GRIDNODEPL_H_
#define OCTI_GRIDGRAPH_GRIDNODEPL_H_

#include "octi/gridgraph/GridEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace gridgraph {

class GridNodePL;
typedef util::graph::Node<GridNodePL, GridEdgePL> GridNode;
typedef util::graph::Edge<GridNodePL, GridEdgePL> GridEdge;

class GridNodePL : util::geograph::GeoNodePL<double> {
 public:
  GridNodePL(){};
  GridNodePL(Point<double> pos);

  const Point<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  GridNode* getParent() const;
  void setParent(GridNode* n);

  GridNode* getPort(size_t i) const;
  void setPort(size_t p, GridNode* n);

  void setXY(size_t x, size_t y);
  size_t getX() const;
  size_t getY() const;

  bool isClosed() const;
  void setClosed(bool c);

 private:
  Point<double> _pos;

  GridNode* _parent;
  GridNode* _ports[8];

  size_t _x, _y;
  bool _closed;
};
}
}

#endif  // OCTI_GRIDGRAPH_GRIDNODEPL_H_
