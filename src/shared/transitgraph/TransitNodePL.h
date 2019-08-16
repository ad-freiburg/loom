// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_TRANSITGRAPH_TRANSITNODEPL_H_
#define SHARED_TRANSITGRAPH_TRANSITNODEPL_H_

#include "shared/transitgraph/TransitEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

#include "transitmap/graph/Node.h"

using util::geo::Point;
using util::graph::Node;
using transitmapper::graph::StationInfo;

namespace shared {
namespace transitgraph {

typedef util::graph::Node<TransitNodePL, TransitEdgePL> TransitNode;
typedef util::graph::Edge<TransitNodePL, TransitEdgePL> TransitEdge;

struct ConnException {
  ConnException(const TransitEdge* from, const TransitEdge* to)
      : fr(from), to(to) {}
  const TransitEdge* fr;
  const TransitEdge* to;
};

class TransitNodePL : util::geograph::GeoNodePL<double> {
 public:
  TransitNodePL(){};
  TransitNodePL(Point<double> pos);

  const Point<double>* getGeom() const;
  void setGeom(const Point<double>& p);
  util::json::Dict getAttrs() const;

  void addStop(StationInfo i);
  const std::vector<StationInfo>& getStops() const;

 private:
  Point<double> _pos;
  std::vector<StationInfo> _is;
};
}
}

#endif  // SHARED_TRANSITGRAPH_TRANSITNODEPL_H_
