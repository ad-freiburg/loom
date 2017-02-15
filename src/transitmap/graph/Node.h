// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_NODE_H_
#define TRANSITMAP_GRAPH_NODE_H_

#include <set>
#include "gtfsparser/gtfs/Stop.h"
#include "./Route.h"
#include "../geo/PolyLine.h"
#include "./OrderingConfiguration.h"
#include "../util/Geo.h"

using namespace gtfsparser;

namespace transitmapper {
namespace graph {

// forward declarations
class Edge;
class Node;
class EdgeTripGeom;
struct RouteOccurance;

// forward declaration of TransitGraph
class TransitGraph;

using util::geo::Point;

struct NodeFront {
  NodeFront(Edge* e, Node* n) : n(n), edge(e) {

  }

  Node* n; // pointer to node here also

  Point getTripOccPos(const Route* r, const Configuration& c) const;
  Point getTripOccPosUnder(const Route* r, const graph::Configuration& c,
    const Edge* e, const std::vector<size_t>* order) const;
  Point getTripPos(const Edge* e, size_t pos, bool inv) const;

  Edge* edge;

  geo::PolyLine geom;
  void setGeom(const geo::PolyLine& g) { geom = g; };

  // TODO
  double refEtgLengthBefExp;
};

struct Partner {
  const NodeFront* front;
  const Edge* edge;
  const EdgeTripGeom* etg;
  const Route* route;
};

struct InnerGeometry {
  InnerGeometry(geo::PolyLine g, const Route* r, const Edge* e)
  : geom(g), route(r), e(e) {};
  geo::PolyLine geom;
  const Route* route;
  const Edge* e;
};

struct StationInfo {
  StationInfo(const std::string& id, const std::string& name) : id(id), name(name) {}
  std::string id, name;
};

class Node {
 public:
  Node(const std::string& id, Point pos);
  Node(const std::string& id, double x, double y);
  Node(const std::string& id, Point pos, StationInfo stop);
  Node(const std::string& id, double x, double y, StationInfo stop);

  ~Node();

  const std::vector<StationInfo>& getStops() const;
  void addStop(StationInfo s);
  const Point& getPos() const;
  void setPos(const Point& p);

  const std::string& getId() const;

  const std::set<Edge*>& getAdjListOut() const {
    return _adjListOut;
  }
  const std::set<Edge*>& getAdjListIn() const {
    return _adjListIn;
  }
  const std::vector<NodeFront>& getMainDirs() const {
    return _mainDirs;
  }
  std::vector<NodeFront>& getMainDirs() {
    return _mainDirs;
  }

  void addMainDir(NodeFront f);

  const NodeFront* getNodeFrontFor(const Edge* e) const;
  double getScore(const graph::Configuration& c) const;
  double getScoreUnder(const graph::Configuration& c, const Edge* e, const graph::Ordering* order) const;
  double getAreaScore(const Configuration& c, const Edge* e, const graph::Ordering* order) const;
  double getAreaScore(const Configuration& c) const;
  std::vector<Partner> getPartners(const NodeFront* f, const RouteOccurance& ro) const;

  std::vector<InnerGeometry> getInnerGeometries(const graph::Configuration& c,
      bool bezier) const;
  std::vector<InnerGeometry> getInnerGeometriesUnder(
      const graph::Configuration& c, bool bezier, const Edge* e,
      const std::vector<size_t>* order) const;

  util::geo::Polygon getConvexFrontHull(double d) const;

  // add edge to this node's adjacency lists
  void addEdge(Edge* e);

  // remove edge from this node's adjacency lists
  void removeEdge(Edge* e);

  double getMaxNodeFrontWidth() const;

 private:
  std::string _id;
  std::set<Edge*> _adjListIn;
  std::set<Edge*> _adjListOut;
  Point _pos;

  std::vector<NodeFront> _mainDirs;

  std::vector<StationInfo> _stops;

  size_t getNodeFrontPos(const NodeFront* a) const;

  geo::PolyLine getInnerBezier(const Configuration& c, const NodeFront& nf,
      const RouteOccurance& tripOcc, const graph::Partner& partner, const Edge* e,
      const std::vector<size_t>* order) const;

  geo::PolyLine getInnerStraightLine(const Configuration& c,
      const NodeFront& nf, const RouteOccurance& tripOcc,
      const graph::Partner& partner, const Edge* e,
      const std::vector<size_t>* order) const;

  friend class TransitGraph;
};
}}

#endif  // TRANSITMAP_GRAPH_NODE_H_
