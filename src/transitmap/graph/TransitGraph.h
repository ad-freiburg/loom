// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_TRANSITGRAPH_H_
#define TRANSITMAP_GRAPH_TRANSITGRAPH_H_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <set>
#include <string>

#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/Route.h"
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/Node.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Penalties.h"
#include "util/geo/Geo.h"

using namespace util;

namespace transitmapper {
namespace graph {

using shared::linegraph::LineNode;
using shared::linegraph::LineEdge;

class TransitGraph : public shared::linegraph::LineGraph {
 public:
  const OrderingConfig& getConfig() const;
  void setConfig(const OrderingConfig&);

  size_t getNumNodes() const;
  size_t getNumNodes(bool topo) const;
  size_t getNumEdges() const;
  size_t getNumRoutes() const;
  size_t getMaxCardinality() const;

  std::vector<shared::linegraph::InnerGeometry> getInnerGeometries(
      const LineNode* n, const graph::OrderingConfig& c, double prec) const;

  Polygon<double> getStationHull(const LineNode* n, double d,
                                 bool simple) const;
  static size_t getConnCardinality(const LineNode* n);

 private:
  mutable double _lastSolveTime;
  mutable size_t _lastSolveTarget;

  OrderingConfig _config;

  shared::linegraph::InnerGeometry getInnerBezier(const LineNode* n,
                                                  const OrderingConfig& cf,
                                                  const Partner& partnerFrom,
                                                  const Partner& partnerTo,
                                                  double prec) const;

  shared::linegraph::InnerGeometry getInnerStraightLine(
      const LineNode* n, const OrderingConfig& c,
      const graph::Partner& partnerFrom, const graph::Partner& partnerTo) const;

  InnerGeometry getTerminusStraightLine(
      const LineNode* n, const OrderingConfig& c,
      const graph::Partner& partnerFrom) const;

  InnerGeometry getTerminusBezier(const LineNode* n, const OrderingConfig& c,
                                  const graph::Partner& partnerFrom,
                                  double prec) const;

  Polygon<double> getConvexFrontHull(const LineNode* n, double d,
                                     bool rectangulize,
                                     bool simpleRenderForTwoEdgeNodes) const;
};
}
}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
