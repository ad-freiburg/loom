// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SKELETONBUILDER_GRAPH_TRANSITGRAPH_H_
#define SKELETONBUILDER_GRAPH_TRANSITGRAPH_H_

#include <proj_api.h>
#include <string>
#include <set>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "transitmap/util/Geo.h"
#include "./Node.h"
#include "./Edge.h"

namespace bg = bgeo;

namespace skeletonbuilder {
namespace graph {

class Graph {

 public:
  explicit Graph(const std::string& name, const std::string& projStr);

  ~Graph();

  Node* addNode(Node* n);
  Edge* addEdge(Node* from, Node* to);
  Edge* getEdge(Node* from, Node* to);

  void deleteNode(Node* n);
  void deleteEdge(Node* from, Node* to);

  const std::set<Node*>& getNodes() const;
  std::set<Node*>* getNodes();

  Node* getNodeByStop(const gtfs::Stop* s) const;
  Node* getNodeByStop(const gtfs::Stop* s, bool getParent) const;

  Node* getNearestNode(const util::geo::Point& p, double maxD) const;

  projPJ getProjection() const;

 private:
  std::string _name;
  std::set<Node*> _nodes;

  projPJ _proj;
};

}}

#endif  // SKELETONBUILDER_GRAPH_TRANSITGRAPH_H_
