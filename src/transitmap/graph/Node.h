// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_NODE_H_
#define TRANSITMAP_GRAPH_NODE_H_

#include <cmath>
#include <set>
#include "shared/linegraph/LineNodePL.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/graph/Route.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"

using namespace util::geo;

namespace transitmapper {
namespace graph {

// forward declarations
class Edge;
class Node;
struct RouteOccurance;

// forward declaration of TransitGraph
class TransitGraph;

struct NodeFront {
  NodeFront(Edge* e, Node* n) : n(n), edge(e) {}

  Node* n;  // pointer to node here also
  DPoint getTripOccPos(const Route* r, const OrderingConfig& c) const;
  DPoint getTripOccPos(const Route* r, const OrderingConfig& c,
                       bool originGeom) const;
  DPoint getTripPos(const Edge* e, size_t pos, bool inv) const;
  DPoint getTripPos(const Edge* e, size_t pos, bool inv, bool originGeom) const;

  double getOutAngle() const;

  Edge* edge;

  // geometry after expansion
  PolyLine<double> geom;

  // geometry before expansion
  PolyLine<double> origGeom;

  void setInitialGeom(const PolyLine<double>& g) {
    geom = g;
    origGeom = g;
  };
  void setGeom(const PolyLine<double>& g) { geom = g; };

  // TODO
  double refEtgLengthBefExp;
};

struct Partner {
  Partner() : front(0), edge(0), route(0){};
  Partner(const NodeFront* f, const Edge* e, const Route* r)
      : front(f), edge(e), route(r){};
  const NodeFront* front;
  const Edge* edge;
  const Route* route;
};

struct InnerGeometry {
  InnerGeometry(PolyLine<double> g, Partner a, Partner b, size_t slotF,
                size_t slotT)
      : geom(g), from(a), to(b), slotFrom(slotF), slotTo(slotT){};
  PolyLine<double> geom;
  Partner from, to;
  size_t slotFrom, slotTo;
};

class Node {
 public:
  Node(const std::string& id, DPoint pos);
  Node(const std::string& id, double x, double y);
  Node(const std::string& id, DPoint pos, shared::linegraph::Station stop);
  Node(const std::string& id, double x, double y,
       shared::linegraph::Station stop);

  ~Node();

  const std::vector<shared::linegraph::Station>& getStops() const;
  void addStop(shared::linegraph::Station s);
  const DPoint& getPos() const;
  void setPos(const DPoint& p);

  const std::string& getId() const;

  const std::set<Edge*>& getAdjListOut() const { return _adjListOut; }
  const std::set<Edge*>& getAdjListIn() const { return _adjListIn; }
  std::set<Edge*> getAdjList() const;

  const std::vector<NodeFront>& getMainDirs() const { return _mainDirs; }
  std::vector<NodeFront>& getMainDirs() { return _mainDirs; }

  void addMainDir(NodeFront f);

  const NodeFront* getNodeFrontFor(const Edge* e) const;

  std::vector<Partner> getPartners(const NodeFront* f,
                                   const RouteOccurance& ro) const;

  std::vector<InnerGeometry> getInnerGeometries(const graph::OrderingConfig& c,
                                                double prec) const;

  size_t getConnCardinality() const;

  Polygon<double> getConvexFrontHull(double d, bool rectangulize,
                                     bool simple) const;

  void generateStationHull(double d, bool useSimple);

  Polygon<double> getStationHull() const;

  // add edge to this node's adjacency lists
  void addEdg(Edge* e);

  // get edge from or to this node, from or to node "other"
  Edge* getEdg(const Node* other) const;

  // remove edge from this node's adjacency lists
  void removeEdge(Edge* e);

  void addRouteConnException(const Route* r, const Edge* edgeA,
                             const Edge* edgeB);
  bool connOccurs(const Route* r, const Edge* edgeA, const Edge* edgeB) const;

  double getMaxNodeFrontWidth() const;
  size_t getMaxNodeFrontCardinality() const;

 private:
  std::string _id;
  std::set<Edge*> _adjListIn;
  std::set<Edge*> _adjListOut;
  DPoint _pos;

  Polygon<double> _stationHull;

  std::vector<NodeFront> _mainDirs;

  std::vector<shared::linegraph::Station> _stops;

  std::map<const Route*, std::map<const Edge*, std::set<const Edge*> > >
      _routeConnExceptions;

  InnerGeometry getInnerBezier(const OrderingConfig& c,
                               const graph::Partner& partnerFrom,
                               const graph::Partner& partnerTo,
                               double prec) const;

  InnerGeometry getInnerStraightLine(const OrderingConfig& c,
                                     const graph::Partner& partnerFrom,
                                     const graph::Partner& partnerTo) const;
  InnerGeometry getTerminusStraightLine(
      const OrderingConfig& c, const graph::Partner& partnerFrom) const;
  InnerGeometry getTerminusBezier(const OrderingConfig& c,
                                  const graph::Partner& partnerFrom,
                                  double prec) const;

  friend class TransitGraph;
};
}
}

#endif  // TRANSITMAP_GRAPH_NODE_H_
