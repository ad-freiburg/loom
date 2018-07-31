// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_NODE_H_
#define TRANSITMAP_GRAPH_NODE_H_

#include <set>
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

  FPoint getTripOccPos(const Route* r, const OrderingConfig& c) const;
  FPoint getTripOccPos(const Route* r, const OrderingConfig& c,
                       bool originGeom) const;
  FPoint getTripPos(const Edge* e, size_t pos, bool inv) const;
  FPoint getTripPos(const Edge* e, size_t pos, bool inv, bool originGeom) const;

  Edge* edge;

  // geometry after expansion
  PolyLine<float> geom;

  // geometry before expansion
  PolyLine<float> origGeom;

  void setInitialGeom(const PolyLine<float>& g) {
    geom = g;
    origGeom = g;
  };
  void setGeom(const PolyLine<float>& g) { geom = g; };

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
  InnerGeometry(PolyLine<float> g, Partner a, Partner b, size_t slotF,
                size_t slotT)
      : geom(g), from(a), to(b), slotFrom(slotF), slotTo(slotT){};
  PolyLine<float> geom;
  Partner from;
  Partner to;
  size_t slotFrom;
  size_t slotTo;
};

struct StationInfo {
  StationInfo(const std::string& id, const std::string& name)
      : id(id), name(name) {}
  std::string id, name;
};

class Node {
 public:
  Node(const std::string& id, FPoint pos);
  Node(const std::string& id, double x, double y);
  Node(const std::string& id, FPoint pos, StationInfo stop);
  Node(const std::string& id, double x, double y, StationInfo stop);

  ~Node();

  const std::vector<StationInfo>& getStops() const;
  void addStop(StationInfo s);
  const FPoint& getPos() const;
  void setPos(const FPoint& p);

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

  Polygon<float> getConvexFrontHull(double d, bool rectangulize,
                                    bool simple) const;

  void generateStationHull(double d, bool useSimple);

  Polygon<float> getStationHull() const;

  // add edge to this node's adjacency lists
  void addEdge(Edge* e);

  // get edge from or to this node, from or to node "other"
  Edge* getEdge(const Node* other) const;

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
  FPoint _pos;

  Polygon<float> _stationHull;

  std::vector<NodeFront> _mainDirs;

  std::vector<StationInfo> _stops;

  std::map<const Route*, std::map<const Edge*, std::set<const Edge*> > >
      _routeConnExceptions;

  size_t getNodeFrontPos(const NodeFront* a) const;

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
