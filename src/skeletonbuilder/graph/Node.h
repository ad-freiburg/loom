// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SKELETONBUILDER_GRAPH_NODE_H_
#define SKELETONBUILDER_GRAPH_NODE_H_

#include <set>
#include "gtfsparser/gtfs/Stop.h"
#include "gtfsparser/gtfs/Route.h"
#include "transitmap/geo/PolyLine.h"
#include "transitmap/util/Geo.h"

using namespace gtfsparser;
using namespace transitmapper;

namespace skeletonbuilder {
namespace graph {

// forward declarations
class Edge;
class Node;
class EdgeTripGeom;

// forward declaration of TransitGraph
class Graph;

using util::geo::Point;

struct OccuringConnection {
  OccuringConnection(const Edge* from, const Edge* to) : from(from), to(to) {}
  const Edge* from;
  const Edge* to;
};

class Node {

 public:
  explicit Node(Point pos);
  Node(double x, double y);
  Node(Point pos, gtfs::Stop* stop);
  Node(double x, double y, gtfs::Stop* stop);

  ~Node();

  const std::set<gtfs::Stop*>& getStops() const;
  void addStop(gtfs::Stop* s);
  const Point& getPos() const;
  void setPos(const Point& p);

  const std::set<Edge*>& getAdjListOut() const;
  const std::set<Edge*>& getAdjListIn() const;

  std::set<Edge*> getAdjList() const;

  // add edge to this node's adjacency lists
  void addEdge(Edge* e);

  // remove edge from this node's adjacency lists
  void removeEdge(Edge* e);

  bool isConnOccuring(const gtfs::Route*, const Edge* from, const Edge* to) const;

  void connOccurs(const gtfs::Route*, const Edge* from, const Edge* to);

  void replaceEdgeInConnections(const Edge* oldE, const Edge* newE);

  void sewConnectionsTogether(Edge* a, Edge* b);

 private:
  std::set<Edge*> _adjListIn;
  std::set<Edge*> _adjListOut;
  Point _pos;

  std::set<gtfs::Stop*> _stops;

  std::map<const gtfs::Route*, std::vector<OccuringConnection> > _occConns;

  friend class Graph;
};
}}

#endif  // SKELETONBUILDER_GRAPH_NODE_H_
