// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_NODEPL_H_
#define OCTI_GRIDGRAPH_NODEPL_H_

#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"
#include "octi/gridgraph/EdgePL.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace gridgraph {

class NodePL : util::geograph::GeoNodePL {
 public:
  NodePL(Point pos);

  const Point* getGeom() const;
  void getAttrs(json::object_t& obj) const;

  Node<NodePL, EdgePL>* getNorthPort() const;
  void setNorthPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getSouthPort() const;
  void setSouthPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getEastPort() const;
  void setEastPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getWestPort() const;
  void setWestPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getSouthEastPort() const;
  void setSouthEastPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getSouthWestPort() const;
  void setSouthWestPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getNorthEastPort() const;
  void setNorthEastPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getNorthWestPort() const;
  void setNorthWestPort(Node<NodePL, EdgePL>* n);

  Node<NodePL, EdgePL>* getPort(size_t i) const;
 private:
  Point _pos;

  Node<NodePL, EdgePL>* _ports[8];
};
}}

#endif  // OCTI_GRIDGRAPH_NODEPL_H_
