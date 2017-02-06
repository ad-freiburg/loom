// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_TRANSITGRAPH_H_
#define TRANSITMAP_GRAPH_TRANSITGRAPH_H_

#include <proj_api.h>
#include <string>
#include <set>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "./../util/Geo.h"
#include "./OrderingConfiguration.h"
#include "./Node.h"
#include "./Edge.h"
#include "./Route.h"

namespace bg = bgeo;

namespace transitmapper {
namespace graph {

class TransitGraph {

 public:
  TransitGraph(const std::string& name, const std::string& projStr);

  ~TransitGraph();

  void addNode(Node* n);
  Edge* addEdge(Node* from, Node* to, geo::PolyLine pl, double w,
    double s);
  Edge* getEdge(Node* from, Node* to);

  void deleteEdge(Node* from, Node* to);

  const std::set<Node*>& getNodes() const;
  std::set<Node*>* getNodes();

  Node* getNodeById(const std::string& id) const;

  Node* getNearestNode(const util::geo::Point& p, double maxD) const;

  projPJ getProjection() const;
  const bgeo::model::box<util::geo::Point>& getBoundingBox() const;

  double getScore() const;
  double getScore(const Configuration& c) const;

  const Configuration& getConfig() const;
  void setConfig(const Configuration&);

  Route* addRoute(const Route* r);
  const Route* getRoute(const std::string& id) const;

  void expandBBox(const Point& p);

 private:
  std::string _name;
  std::set<Node*> _nodes;
  std::map<std::string, const Route*> _routes;

  Configuration _config;

  bgeo::model::box<util::geo::Point> _bbox;

  projPJ _proj;
};

}}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
