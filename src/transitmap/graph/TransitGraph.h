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

class TransitGraph {
 public:
  TransitGraph();
  ~TransitGraph();

  // +
  void addNd(Node* n);

  // +
  Edge* addEdg(Node* from, Node* to, geo::PolyLine<double> pl, double w,
               double s);
  // +
  Edge* getEdg(Node* from, Node* to);

  // +
  void delEdg(Node* from, Node* to);

  // +
  const std::set<Node*>& getNds() const;
  // +
  std::set<Node*>* getNds();

  // +
  void addRoute(const shared::linegraph::Route* r);

  Node* getNodeById(const std::string& id) const;

  // +
  const geo::DBox& getBBox() const;

  // +
  size_t maxDeg() const;

  const OrderingConfig& getConfig() const;
  void setConfig(const OrderingConfig&);

  // +
  const shared::linegraph::Route* getRoute(const std::string& id) const;

  // +
  void expandBBox(const geo::DPoint& p);

  size_t getNumNodes() const;
  size_t getNumNodes(bool topo) const;
  size_t getNumEdges() const;
  size_t getNumRoutes() const;
  size_t getMaxCardinality() const;

  // +
  bool readFromJson(std::istream* s);

  std::vector<InnerGeometry> getInnerGeometries(const Node* n,
                                                const graph::OrderingConfig& c,
                                                double prec) const;

  Polygon<double> getStationHull(const Node* n, double d, bool simple) const;

 private:
  std::set<Node*> _nodes;
  std::map<std::string, const shared::linegraph::Route*> _routes;

  mutable double _lastSolveTime;
  mutable size_t _lastSolveTarget;

  OrderingConfig _config;

  geo::DBox _bbox;

  InnerGeometry getInnerBezier(const Node* n, const OrderingConfig& cf,
                               const Partner& partnerFrom,
                               const Partner& partnerTo, double prec) const;

  InnerGeometry getInnerStraightLine(const Node* n, const OrderingConfig& c,
                                     const graph::Partner& partnerFrom,
                                     const graph::Partner& partnerTo) const;

  InnerGeometry getTerminusStraightLine(
      const Node* n, const OrderingConfig& c,
      const graph::Partner& partnerFrom) const;

  InnerGeometry getTerminusBezier(const Node* n, const OrderingConfig& c,
                                  const graph::Partner& partnerFrom,
                                  double prec) const;

  Polygon<double> getConvexFrontHull(const Node* n, double d, bool rectangulize,
                                     bool simpleRenderForTwoEdgeNodes) const;
};
}
}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
