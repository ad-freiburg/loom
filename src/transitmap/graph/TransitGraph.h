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
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Penalties.h"
#include "util/geo/Geo.h"

using namespace util;

namespace transitmapper {
namespace graph {

class TransitGraph : public shared::linegraph::LineGraph {
 public:
  TransitGraph(double defaultLineWidth, double defaultLineSpacing)
      : _defWidth(defaultLineWidth), _defSpacing(defaultLineSpacing){};
  const OrderingConfig& getConfig() const;
  void setConfig(const OrderingConfig&);

  size_t getNumNodes() const;
  size_t getNumNodes(bool topo) const;
  size_t getNumEdges() const;
  size_t getNumRoutes() const;
  size_t getMaxCardinality() const;

  std::vector<shared::linegraph::InnerGeometry> getInnerGeometries(
      const shared::linegraph::LineNode* n, const graph::OrderingConfig& c,
      double prec) const;

  util::geo::Polygon<double> getStationHull(
      const shared::linegraph::LineNode* n, double d, bool simple) const;
  static size_t getConnCardinality(const shared::linegraph::LineNode* n);
  double getTotalWidth(const shared::linegraph::LineEdge* e) const;
  double getWidth(const shared::linegraph::LineEdge* e) const;
  double getSpacing(const shared::linegraph::LineEdge* e) const;

  double getMaxNodeFrontWidth(const shared::linegraph::LineNode* n) const;
  size_t getMaxNodeFrontCard(const shared::linegraph::LineNode* n) const;

  double getNodeFrontAngle(const shared::linegraph::NodeFront& nf) const;

  util::geo::DPoint getTripOccPos(const shared::linegraph::NodeFront& nf,
                                  const shared::linegraph::Route* r,
                                  const OrderingConfig& c, bool origGeom) const;

  util::geo::DPoint getTripPos(const shared::linegraph::NodeFront& nf,
                               const shared::linegraph::LineEdge* e, size_t pos,
                               bool inv, bool origG) const;

 private:
  mutable double _lastSolveTime;
  mutable size_t _lastSolveTarget;

  double _defWidth, _defSpacing;

  OrderingConfig _config;

  shared::linegraph::InnerGeometry getInnerBezier(
      const shared::linegraph::LineNode* n, const OrderingConfig& cf,
      const shared::linegraph::Partner& partnerFrom,
      const shared::linegraph::Partner& partnerTo, double prec) const;

  shared::linegraph::InnerGeometry getInnerStraightLine(
      const shared::linegraph::LineNode* n, const OrderingConfig& c,
      const shared::linegraph::Partner& partnerFrom,
      const shared::linegraph::Partner& partnerTo) const;

  shared::linegraph::InnerGeometry getTerminusStraightLine(
      const shared::linegraph::LineNode* n, const OrderingConfig& c,
      const shared::linegraph::Partner& partnerFrom) const;

  shared::linegraph::InnerGeometry getTerminusBezier(
      const shared::linegraph::LineNode* n, const OrderingConfig& c,
      const shared::linegraph::Partner& partnerFrom, double prec) const;

  util::geo::Polygon<double> getConvexFrontHull(
      const shared::linegraph::LineNode* n, double d, bool rectangulize,
      bool simpleRenderForTwoEdgeNodes) const;
};
}
}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
