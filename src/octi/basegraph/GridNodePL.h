// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_BASEGRAPH_GRIDNODEPL_H_
#define OCTI_BASEGRAPH_GRIDNODEPL_H_

#include "octi/basegraph/GridEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace basegraph {

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

  bool isSink() const;
  void setSink();

  bool isSettled() const;
  void setSettled(bool c);

  void setId(size_t id);
  size_t getId() const;

 private:
  Point<double> _pos;

  GridNode* _parent;
  GridNode* _ports[8];

  uint32_t _x, _y;
  uint32_t _id;
  bool _closed : 1;
  bool _sink : 1;
  bool _settled : 1;
};
}  // namespace basegraph
}  // namespace octi

#endif  // OCTI_BASEGRAPH_GRIDNODEPL_H_
