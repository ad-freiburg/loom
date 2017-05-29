// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_TRANSITGRAPH_H_
#define TRANSITMAP_GRAPH_TRANSITGRAPH_H_

#include <proj_api.h>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <set>
#include <string>

#include "./Edge.h"
#include "./Node.h"
#include "./OrderingConfiguration.h"
#include "./Route.h"
#include "pbutil/geo/Geo.h"

namespace bg = bgeo;
using namespace pbutil::geo;

namespace transitmapper {
namespace graph {

class TransitGraph {
 public:
  TransitGraph(const std::string& name, const std::string& projStr);

  ~TransitGraph();

  void addNode(Node* n);
  Edge* addEdge(Node* from, Node* to, PolyLine pl, double w, double s);
  Edge* getEdge(Node* from, Node* to);

  void deleteEdge(Node* from, Node* to);

  const std::set<Node*>& getNodes() const;
  std::set<Node*>* getNodes();

  Node* getNodeById(const std::string& id) const;

  Node* getNearestNode(const Point& p, double maxD) const;

  projPJ getProjection() const;
  const bgeo::model::box<Point>& getBoundingBox() const;

  double getScore(double inStatPen, double sameSegCrossPen,
                  double diffSegCrossPen, double splitPen) const;
  double getScore(double inStatPen, double sameSegCrossPen,
                  double diffSegCrossPen, double splitPen,
                  const Configuration& c) const;

  double getCrossScore(double inStatPen, double sameSegCrossPen,
                       double diffSegCrossPen) const;
  double getCrossScore(double inStatPen, double sameSegCrossPen,
                       double diffSegCrossPen, const Configuration& c) const;

  double getSeparationScore(double inStatPen, double splitPen) const;
  double getSeparationScore(double inStatPen, double splitPen,
                            const Configuration& c) const;

  size_t getNumCrossings() const;
  size_t getNumCrossings(const Configuration& c) const;

  size_t getNumSeparations() const;
  size_t getNumSeparations(const Configuration& c) const;

  double getNumPossSolutions() const;

  const Configuration& getConfig() const;
  void setConfig(const Configuration&);

  Route* addRoute(const Route* r);
  const Route* getRoute(const std::string& id) const;

  void expandBBox(const Point& p);

  size_t getNumNodes() const;
  size_t getNumNodes(bool topo) const;
  size_t getNumEdges() const;
  size_t getNumRoutes() const;
  size_t getMaxCardinality() const;

  const std::string& getName() const;

 private:
  std::string _name;
  std::set<Node*> _nodes;
  std::map<std::string, const Route*> _routes;

  Configuration _config;

  bgeo::model::box<Point> _bbox;

  projPJ _proj;
};
}
}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
