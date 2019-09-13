// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_TRANSITGRAPH_TRANSITGRAPH_H_
#define SHARED_TRANSITGRAPH_TRANSITGRAPH_H_

#include "shared/transitgraph/TransitEdgePL.h"
#include "shared/transitgraph/TransitNodePL.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/UndirGraph.h"

using util::geo::Grid;
using util::geo::Point;

namespace shared {
namespace transitgraph {

typedef util::graph::Node<TransitNodePL, TransitEdgePL> TransitNode;
typedef util::graph::Edge<TransitNodePL, TransitEdgePL> TransitEdge;

typedef std::pair<TransitEdge*, TransitEdge*> TransitEdgePair;

typedef Grid<TransitNode*, util::geo::Point, double> NodeGrid;
typedef Grid<TransitEdge*, util::geo::Line, double> EdgeGrid;

struct ISect {
  TransitEdge *a, *b;
  util::geo::LinePoint<double> bp;
};

class TransitGraph
    : public util::graph::UndirGraph<TransitNodePL, TransitEdgePL> {
 public:
  TransitGraph();

  void readFromJson(std::istream* s);
  void readFromDot(std::istream* s);

  const util::geo::Box<double>& getBBox() const;
  void topologizeIsects();

  static TransitNode* sharedNode(const TransitEdge* a, const TransitEdge* b);
  static std::vector<RouteOcc> getCtdRoutesIn(const Route* r,
                                              const TransitNode* dir,
                                              const TransitEdge* fromEdge,
                                              const TransitEdge* toEdge);

  static std::vector<RouteOcc> getCtdRoutesIn(const TransitEdge* fromEdge,
                                       const TransitEdge* toEdge);

 private:
  util::geo::Box<double> _bbox;

  ISect getNextIntersection();

  void buildGrids();

  void addRoute(const Route* r);
  const Route* getRoute(const std::string& id) const;

  void expandBBox(const util::geo::Point<double>& p);

  std::set<TransitEdge*> proced;
  std::map<std::string, const Route*> _routes;

  NodeGrid _nodeGrid;
  EdgeGrid _edgeGrid;
};

}  // transitgraph
}  // shared

#endif  // SHARED_TRANSITGRAPH_TRANSITGRAPH_H_
