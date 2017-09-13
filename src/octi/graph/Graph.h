// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_GRAPH_H_
#define OCTI_GRAPH_GRAPH_H_

#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/Graph.h"
#include "octi/graph/NodePL.h"
#include "octi/graph/EdgePL.h"

using util::geo::Grid;
using util::geo::Point;

namespace octi {
namespace graph {

struct ISect {
  util::graph::Edge<NodePL, EdgePL>* a;
  util::graph::Edge<NodePL, EdgePL>* b;

  util::geo::LinePoint bp;
};

class Graph : public util::graph::Graph<NodePL, EdgePL> {
 public:
  Graph(std::istream* s);
  Graph();

  const util::geo::Box& getBBox() const;

 private:
  util::geo::Box _bbox;

  void topologizeIsects();
  ISect getNextIntersection();

  void addRoute(const Route* r);
  const Route* getRoute(const std::string& id) const;

  void expandBBox(const Point& p);

  std::set<util::graph::Edge<NodePL, EdgePL>*> proced;
  std::map<std::string, const Route*> _routes;
};

}  // graph
}  // octi

#endif  // OCTI_GRAPH_GRAPH_H_
