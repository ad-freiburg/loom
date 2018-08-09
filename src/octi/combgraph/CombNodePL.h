// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_COMBNODEPL_H_
#define OCTI_COMBGRAPH_COMBNODEPL_H_

#include "octi/combgraph/CombEdgePL.h"
#include "octi/transitgraph/EdgeOrdering.h"
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
  CombNodePL(octi::transitgraph::TransitNode* parent);

  const Point<double>* getGeom() const;
  util::json::Dict getAttrs() const;
  octi::transitgraph::TransitNode* getParent() const;

  const transitgraph::EdgeOrdering& getEdgeOrdering();
  void setEdgeOrdering(const transitgraph::EdgeOrdering& e);
  size_t getRouteNumber() const;
  void setRouteNumber(size_t n);
  std::string toString() const;

 private:
  octi::transitgraph::TransitNode* _parent;
  size_t _routeNumber;
  transitgraph::EdgeOrdering _ordering;
};
}
}

#endif  // OCTI_COMBGRAPH_COMBNODEPL_H_
