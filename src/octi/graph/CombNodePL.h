// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_COMBNODEPL_H_
#define OCTI_GRAPH_COMBNODEPL_H_

#include "octi/graph/CombEdgePL.h"
#include "octi/graph/EdgeOrdering.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace graph {

typedef util::graph::Node<NodePL, EdgePL> Node;

// forward declaration
class CombNodePL;

class CombNodePL : util::geograph::GeoNodePL {
 public:
  CombNodePL(octi::graph::Node* parent);

  const Point* getGeom() const;
  void getAttrs(json::object_t& obj) const;
  Node* getParent() const;

  const EdgeOrdering& getEdgeOrdering();
  void setEdgeOrdering(const EdgeOrdering& e);
  size_t getRouteNumber() const;
  void setRouteNumber(size_t n);
  std::string toString() const;

 private:
  octi::graph::Node* _parent;
  size_t _routeNumber;
  EdgeOrdering _ordering;

};
}
}

#endif  // OCTI_GRAPH_COMBNODEPL_H_
