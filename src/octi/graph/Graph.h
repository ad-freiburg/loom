// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_GRAPH_H_
#define OCTI_GRAPH_GRAPH_H_

#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/UndirGraph.h"
#include "octi/graph/NodePL.h"
#include "octi/graph/EdgePL.h"

using util::geo::Grid;
using util::geo::Point;

namespace octi {
namespace graph {

typedef Grid<util::graph::Node<NodePL, EdgePL>*, Point, double> NodeGrid;
typedef Grid<util::graph::Edge<NodePL, EdgePL>*, Line, double> EdgeGrid;

struct ISect {
  util::graph::Edge<NodePL, EdgePL>* a;
  util::graph::Edge<NodePL, EdgePL>* b;

  util::geo::LinePoint<double> bp;
};

class Graph : public util::graph::UndirGraph<NodePL, EdgePL> {
 public:
  Graph();

  void readFromJson(std::istream* s);
  void readFromDot(std::istream* s);

  const util::geo::Box<double>& getBBox() const;
  void topologizeIsects();

 private:
  util::geo::Box<double> _bbox;

  ISect getNextIntersection();

  void buildGrids();

  void addRoute(const Route* r);
  const Route* getRoute(const std::string& id) const;

  void expandBBox(const Point<double>& p);

  std::set<util::graph::Edge<NodePL, EdgePL>*> proced;
  std::map<std::string, const Route*> _routes;

  NodeGrid _nodeGrid;
  EdgeGrid _edgeGrid;
};

}  // graph
}  // octi

#endif  // OCTI_GRAPH_GRAPH_H_
