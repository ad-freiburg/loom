// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_NODE_H_
#define TRANSITMAP_GRAPH_NODE_H_

#include <set>
#include "gtfsparser/gtfs/Stop.h"
#include "gtfsparser/gtfs/Route.h"
#include "../geo/PolyLine.h"
#include "../util/Geo.h"

using namespace gtfsparser;

namespace transitmapper {
namespace graph {

// forward declaration of Edge
class Edge;

class Node;

// forward declaration of EdgeTripGeometry
class EdgeTripGeom;

// forward declaration of TransitGraph
class TransitGraph;

struct NodeFront {
  NodeFront(double angle, Edge* e, Node* n) : angle(angle), n(n) {
    addEdge(e);
  }
  double angle;

  Node* n; // pointer to node here also

  std::vector<Edge*> edges;

  util::geo::Point getTripOccPos(const gtfs::Route*) const;

  geo::PolyLine geom;
  void setGeom(const geo::PolyLine& g) { geom = g; };
  void addEdge(Edge* e) { edges.push_back(e); }
};

struct Partner {
  const NodeFront* front;
  const Edge* edge;
  const EdgeTripGeom* etg;
  const gtfs::Route* route;
};

struct InnerGeometry {
  InnerGeometry(geo::PolyLine g, const gtfs::Route* r, const EdgeTripGeom* etg)
  : geom(g), route(r), etg(etg) {};
  geo::PolyLine geom;
  const gtfs::Route* route;
  const EdgeTripGeom* etg;
};

class Node {

 public:
  explicit Node(util::geo::Point pos);
  Node(double x, double y);
  Node(util::geo::Point pos, gtfs::Stop* stop);
  Node(double x, double y, gtfs::Stop* stop);

  ~Node();

  const std::set<gtfs::Stop*>& getStops() const;
  void addStop(gtfs::Stop* s);
  const util::geo::Point& getPos() const;
  void setPos(const util::geo::Point& p);

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
  double getScore() const;
  std::vector<Partner> getPartner(const NodeFront* f, const gtfs::Route* r) const;

  std::vector<InnerGeometry> getInnerGeometries() const;

 protected:
  // add edge to this node's adjacency lists
  void addEdge(Edge* e);

  // remove edge from this node's adjacency lists
  void removeEdge(Edge* e);

 private:
  std::set<Edge*> _adjListIn;
  std::set<Edge*> _adjListOut;
  util::geo::Point _pos;

  std::vector<NodeFront> _mainDirs;

  std::set<gtfs::Stop*> _stops;

  friend class TransitGraph;
};
}}

#endif  // TRANSITMAP_GRAPH_NODE_H_
