// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_TRANSITGRAPH_H_
#define TRANSITMAP_GRAPH_TRANSITGRAPH_H_

#include <string>
#include <set>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "node.h"
#include "edge.h"

namespace bg = boost::geometry;

namespace transitmapper {
namespace graph {

class TransitGraph {

 public:
  explicit TransitGraph(const std::string& name);

  ~TransitGraph(); // TODO: free memory

  Node* addNode(Node* n);
  Edge* addEdge(Node* from, Node* to);

  void deleteNode(Node* n);
  bool containsNode(Node* n) const;

  const std::set<Node*>& getNodes() const;
  std::set<Node*>* getNodes();

 private:
  std::string _name;
  std::set<Node*> _nodes;
};

}}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
