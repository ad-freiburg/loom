// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_NODEPL_H_
#define OCTI_GRAPH_NODEPL_H_

#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"
#include "octi/graph/EdgePL.h"

#include "transitmap/graph/Node.h"

using util::geo::Point;
using util::graph::Node;
using transitmapper::graph::StationInfo;

namespace octi {
namespace graph {

class NodePL : util::geograph::GeoNodePL {
 public:
  NodePL(Point pos);

  const Point* getGeom() const;
  void setGeom(const Point& p);
  void getAttrs(json::object_t& obj) const;

  void addStop(StationInfo i);

 private:
  Point _pos;
  std::vector<StationInfo> _is;
};
}}

#endif  // OCTI_GRAPH_NODEPL_H_
