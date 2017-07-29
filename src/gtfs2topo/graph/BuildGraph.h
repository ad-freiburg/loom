// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_GRAPH_BUILDGRAPH_H_
#define GTFS2TOPO_GRAPH_BUILDGRAPH_H_

#include <proj_api.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <set>
#include <string>

#include "gtfs2topo/graph/Edge.h"
#include "gtfs2topo/graph/Node.h"
#include "util/geo/Geo.h"

namespace gtfs2topo {
namespace graph {

class BuildGraph {
 public:
  explicit BuildGraph(const std::string& name, const std::string& projStr);
  ~BuildGraph();

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
}
}

#endif  // GTFS2TOPO_GRAPH_BUILDGRAPH_H_
