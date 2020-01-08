// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINEGRAPH_H_
#define SHARED_LINEGRAPH_LINEGRAPH_H_

#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineNodePL.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/UndirGraph.h"

using util::geo::Grid;
using util::geo::Point;

namespace shared {
namespace linegraph {

typedef util::graph::Node<LineNodePL, LineEdgePL> LineNode;
typedef util::graph::Edge<LineNodePL, LineEdgePL> LineEdge;

typedef std::pair<LineEdge*, LineEdge*> LineEdgePair;

typedef Grid<LineNode*, util::geo::Point, double> NodeGrid;
typedef Grid<LineEdge*, util::geo::Line, double> EdgeGrid;

struct ISect {
  LineEdge *a, *b;
  util::geo::LinePoint<double> bp;
};

class LineGraph
    : public util::graph::UndirGraph<LineNodePL, LineEdgePL> {
 public:
  LineGraph();

  void readFromJson(std::istream* s);
  void readFromDot(std::istream* s);

  const util::geo::Box<double>& getBBox() const;
  void topologizeIsects();

  size_t maxDeg() const;

  // TODO: make the following functions private
  void addRoute(const Route* r);
  const Route* getRoute(const std::string& id) const;
  void expandBBox(const util::geo::Point<double>& p);
  //

  static LineNode* sharedNode(const LineEdge* a, const LineEdge* b);
  static std::vector<RouteOcc> getCtdRoutesIn(const Route* r,
                                              const LineNode* dir,
                                              const LineEdge* fromEdge,
                                              const LineEdge* toEdge);

  static std::vector<RouteOcc> getCtdRoutesIn(const LineEdge* fromEdge,
                                              const LineEdge* toEdge);
  static size_t getLDeg(const LineNode* nd);
  static size_t getMaxLineNum(const LineNode* nd);


 private:
  util::geo::Box<double> _bbox;

  ISect getNextIntersection();

  void buildGrids();

  std::set<LineEdge*> proced;
  std::map<std::string, const Route*> _routes;

  NodeGrid _nodeGrid;
  EdgeGrid _edgeGrid;
};

}  // linegraph
}  // shared

#endif  // SHARED_LINEGRAPH_LINEGRAPH_H_
