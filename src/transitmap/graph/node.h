// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_NODE_H_
#define TRANSITMAP_GRAPH_NODE_H_

#include <set>
#include "gtfsparser/gtfs/stop.h"
#include "../util/Geo.h"

using namespace gtfsparser;

namespace transitmapper {
namespace graph {

// forward declaration of Edge
class Edge;

// forward declaration of TransitGraph
class TransitGraph;

struct NodeFront {
  NodeFront(double angle, Edge* e) : angle(angle) {
    addEdge(e);
  }
  double angle;
  std::vector<Edge*> edges;

  void addEdge(Edge* e) { edges.push_back(e); }
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

  void addMainDir(NodeFront f);

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
