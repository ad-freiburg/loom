// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_COMBNODEPL_H_
#define OCTI_COMBGRAPH_COMBNODEPL_H_

#include "octi/combgraph/CombEdgePL.h"
#include "octi/combgraph/EdgeOrdering.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace combgraph {

class CombNodePL : util::geograph::GeoNodePL<double> {
 public:
  CombNodePL(){};
  CombNodePL(shared::linegraph::LineNode* parent);

  const Point<double>* getGeom() const;
  util::json::Dict getAttrs() const;
  shared::linegraph::LineNode* getParent() const;

  const combgraph::EdgeOrdering& getEdgeOrdering();
  void setEdgeOrdering(const combgraph::EdgeOrdering& e);
  size_t getRouteNumber() const;
  void setRouteNumber(size_t n);
  std::string toString() const;

 private:
  shared::linegraph::LineNode* _parent;
  size_t _routeNumber;
  combgraph::EdgeOrdering _ordering;
};
}
}

#endif  // OCTI_COMBGRAPH_COMBNODEPL_H_
