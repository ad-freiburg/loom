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

class Node {

 public:
  explicit Node(util::geo::Point pos);
  Node(double x, double y);
  Node(util::geo::Point pos, gtfs::Stop* stop);
  Node(double x, double y, gtfs::Stop* stop);

  virtual ~Node();

  virtual gtfs::Stop* getStop() const;
  const util::geo::Point& getPos() const;

 protected:
  // add edge to this node's adjacency lists
  void addEdge(Edge* e);

  // remove edge from this node's adjacency lists
  void removeEdge(Edge* e);

 private:
  std::set<Edge*> _adjListIn;
  std::set<Edge*> _adjListOut;
  util::geo::Point _pos;

  gtfs::Stop* _stop;

  friend class TransitGraph;
};
}}

#endif  // TRANSITMAP_GRAPH_NODE_H_
