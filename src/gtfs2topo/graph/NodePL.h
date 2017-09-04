// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_GRAPH_NODE_H_
#define GTFS2TOPO_GRAPH_NODE_H_

#include <set>
#include "ad/cppgtfs/gtfs/Stop.h"
#include "ad/cppgtfs/gtfs/Route.h"
#include "util/geo/PolyLine.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "gtfs2topo/graph/BuildGraph.h"

using namespace ad::cppgtfs;

namespace gtfs2topo {
namespace graph {


using util::geo::Point;

struct OccuringConnection {
  OccuringConnection(const Edge* from, const Edge* to) : from(from), to(to) {}
  const Edge* from;
  const Edge* to;
};

class NodePL : util::geograph::GeoNodePL {
 public:
  NodePL(Point pos);
  NodePL(double x, double y);
  NodePL(Point pos, const gtfs::Stop* stop);
  NodePL(double x, double y, const gtfs::Stop* stop);

  const std::set<const gtfs::Stop*>& getStops() const;
  void addStop(const gtfs::Stop* s);
  const Point& getPos() const;
  void setPos(const Point& p);

  bool isConnOccuring(const gtfs::Route*, const Edge* from, const Edge* to) const;

  void connOccurs(const gtfs::Route*, const Edge* from, const Edge* to);

  void replaceEdgeInConnections(const Edge* oldE, const Edge* newE);

  void sewConnectionsTogether(Edge* a, Edge* b);

  std::vector<const Edge*> getConnectingEdgesFor(const gtfs::Route* to, Edge* a) const;

  const std::map<const gtfs::Route*, std::vector<OccuringConnection> >& getOccuringConnections() const;

  const util::geo::Point* getGeom() const;
  void getAttrs(json::object_t& obj) const;

  void setNode(const Node* n);

 private:
  Point _pos;
  const Node* _n; // backpointer to node

  std::set<const gtfs::Stop*> _stops;
  std::map<const gtfs::Route*, std::vector<OccuringConnection> > _occConns;
};
}}

#endif  // GTFS2TOPO_GRAPH_NODE_H_
