// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_RENDERGRAPH_H_
#define TRANSITMAP_GRAPH_RENDERGRAPH_H_

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <set>
#include <string>

#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/Line.h"
#include "transitmap/graph/OrderCfg.h"
#include "transitmap/graph/Penalties.h"
#include "util/geo/Geo.h"

namespace transitmapper {
namespace graph {

class RenderGraph : public shared::linegraph::LineGraph {
 public:
  RenderGraph(double defLineWidth, double defLineSpace)
      : _defWidth(defLineWidth), _defSpacing(defLineSpace){};

  const OrderCfg& getConfig() const;
  void setConfig(const OrderCfg&);

  size_t numEdgs() const;

  std::vector<shared::linegraph::InnerGeom> innerGeoms(
      const shared::linegraph::LineNode* n, const OrderCfg& c,
      double prec) const;

  util::geo::Polygon<double> getStationHull(
      const shared::linegraph::LineNode* n, double d, bool simple) const;

  // TODO: maybe move this to LineGraph?
  static size_t getConnCardinality(const shared::linegraph::LineNode* n);

  double getTotalWidth(const shared::linegraph::LineEdge* e) const;

  double getWidth(const shared::linegraph::LineEdge* e) const;
  double getSpacing(const shared::linegraph::LineEdge* e) const;

  double getMaxNdFrontWidth(const shared::linegraph::LineNode* n) const;
  size_t getMaxNdFrontCard(const shared::linegraph::LineNode* n) const;

  util::geo::DPoint linePosOn(const shared::linegraph::NodeFront& nf,
                              const shared::linegraph::Line* r,
                              const OrderCfg& c, bool origGeom) const;

  util::geo::DPoint linePosOn(const shared::linegraph::NodeFront& nf,
                              const shared::linegraph::LineEdge* e, size_t pos,
                              bool inv, bool origG) const;

 private:
  double _defWidth, _defSpacing;

  OrderCfg _config;

  shared::linegraph::InnerGeom getInnerBezier(
      const shared::linegraph::LineNode* n, const OrderCfg& cf,
      const shared::linegraph::Partner& partnerFrom,
      const shared::linegraph::Partner& partnerTo, double prec) const;

  shared::linegraph::InnerGeom getInnerStraightLine(
      const shared::linegraph::LineNode* n, const OrderCfg& c,
      const shared::linegraph::Partner& partnerFrom,
      const shared::linegraph::Partner& partnerTo) const;

  shared::linegraph::InnerGeom getTerminusStraightLine(
      const shared::linegraph::LineNode* n, const OrderCfg& c,
      const shared::linegraph::Partner& partnerFrom) const;

  shared::linegraph::InnerGeom getTerminusBezier(
      const shared::linegraph::LineNode* n, const OrderCfg& c,
      const shared::linegraph::Partner& partnerFrom, double prec) const;

  util::geo::Polygon<double> getConvexFrontHull(
      const shared::linegraph::LineNode* n, double d, bool rectangulize,
      bool simpleRenderForTwoEdgeNodes) const;
};
}
}

#endif  // TRANSITMAP_GRAPH_RENDERGRAPH_H_
