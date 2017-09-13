// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_COMBNODEPL_H_
#define OCTI_GRAPH_COMBNODEPL_H_

#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"
#include "octi/graph/CombEdgePL.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace graph {

typedef util::graph::Node<NodePL, EdgePL> Node;

class CombNodePL : util::geograph::GeoNodePL {
 public:
  CombNodePL(octi::graph::Node* parent);

  const Point* getGeom() const;
  void getAttrs(json::object_t& obj) const;
  Node* getParent() const;

 private:
  octi::graph::Node* _parent;
};
}}

#endif  // OCTI_GRAPH_COMBNODEPL_H_
