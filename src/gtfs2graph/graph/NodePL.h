// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2GRAPH_GRAPH_NODE_H_
#define GTFS2GRAPH_GRAPH_NODE_H_

#include <set>
#include "ad/cppgtfs/gtfs/Stop.h"
#include "ad/cppgtfs/gtfs/Route.h"
#include "util/geo/PolyLine.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "gtfs2graph/graph/BuildGraph.h"

using namespace ad::cppgtfs;
using util::geo::DPoint;

namespace gtfs2graph {
namespace graph {

struct OccuringConnection {
  OccuringConnection(const Edge* from, const Edge* to) : from(from), to(to) {}
  const Edge* from;
  const Edge* to;
};

class NodePL : util::geograph::GeoNodePL<double> {
 public:
  NodePL() {};
  NodePL(DPoint pos);
  NodePL(double x, double y);
  NodePL(DPoint pos, const gtfs::Stop* stop);
  NodePL(double x, double y, const gtfs::Stop* stop);

  const std::set<const gtfs::Stop*>& getStops() const;
  void addStop(const gtfs::Stop* s);
  const DPoint& getPos() const;
  void setPos(const DPoint& p);

  bool isConnOccuring(const gtfs::Route*, const Edge* from, const Edge* to) const;

  void connOccurs(const gtfs::Route*, const Edge* from, const Edge* to);

  std::vector<const Edge*> getConnectingEdgesFor(const gtfs::Route* to, Edge* a) const;

  const std::map<const gtfs::Route*, std::vector<OccuringConnection> >& getOccuringConnections() const;

  const util::geo::DPoint* getGeom() const;
  util::json::Dict getAttrs() const;

  void setNode(const Node* n);

 private:
  DPoint _pos;
  const Node* _n; // backpointer to node

  std::set<const gtfs::Stop*> _stops;
  std::map<const gtfs::Route*, std::vector<OccuringConnection> > _occConns;
};
}}

#endif  // GTFS2GRAPH_GRAPH_NODE_H_
