// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_TRANSITGRAPH_H_
#define TRANSITMAP_GRAPH_TRANSITGRAPH_H_

#include <proj_api.h>
#include <string>
#include <set>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "../util/Geo.h"
#include "node.h"
#include "edge.h"

namespace bg = boost::geometry;

namespace transitmapper {
namespace graph {

class TransitGraph {

 public:
  explicit TransitGraph(const std::string& name, const std::string& projStr);

  ~TransitGraph(); // TODO: free memory

  Node* addNode(Node* n);
  Edge* addEdge(Node* from, Node* to);
  Edge* getEdge(Node* from, Node* to);

  void deleteNode(Node* n);
  bool containsNode(Node* n) const;

  const std::set<Node*>& getNodes() const;
  std::set<Node*>* getNodes();

  Node* getNodeByStop(const gtfs::Stop* s) const;
  Node* getNodeByStop(const gtfs::Stop* s, bool getParent) const;

  projPJ getProjection() const;
  const boost::geometry::model::box<util::geo::Point>& getBoundingBox() const;
 private:
  std::string _name;
  std::set<Node*> _nodes;

  boost::geometry::model::box<util::geo::Point> _bbox;

  projPJ _proj;
};

}}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
