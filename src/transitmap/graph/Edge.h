// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_EDGE_H_
#define TRANSITMAP_GRAPH_EDGE_H_

#include <vector>
#include "Node.h"
#include "shared/linegraph/Route.h"
#include "transitmap/style/LineStyle.h"
#include "util/Nullable.h"
#include "util/geo/PolyLine.h"

using std::exception;
using std::string;
using namespace util::geo;

namespace transitmapper {
namespace graph {

using style::LineStyle;

// forward declaration of Node
class Node;

struct RouteOccurance {
  RouteOccurance(const shared::linegraph::Route* r, const Node* dir)
      : route(r), direction(dir) {}
  RouteOccurance(const shared::linegraph::Route* r, const Node* dir,
                 const style::LineStyle& ls)
      : route(r), direction(dir), style(ls) {}
  const shared::linegraph::Route* route;
  const Node* direction;  // 0 if in both directions

  util::Nullable<style::LineStyle> style;

  bool operator==(const RouteOccurance& b) const { return b.route == route; }
};

typedef std::pair<RouteOccurance*, size_t> RouteOccWithPos;

class Edge {
 public:
  Edge(Node* from, Node* to, PolyLine<double> pl, double w, double s);

  Node* getFrom() const;
  Node* getTo() const;

  Node* getOther(const Node* n) const;

  void setFrom(Node* from);
  void setTo(Node* to);

  void addRoute(const shared::linegraph::Route* r, const Node* dir,
                const LineStyle& ls);
  void addRoute(const shared::linegraph::Route* r, const Node* dir);

  const PolyLine<double>& getGeom() const;
  void setGeom(const PolyLine<double>& p);

  const std::vector<RouteOccurance>& getRoutes() const;
  std::vector<RouteOccurance>* getRoutes();

  size_t getRoutePosUnder(const shared::linegraph::Route* r,
                          const std::vector<size_t> ordering) const;

  size_t getRoutePos(const shared::linegraph::Route* r) const;

  RouteOccurance* getRoute(const shared::linegraph::Route* r) const;

  std::vector<RouteOccurance> getCtdRoutesIn(const Node* n,
                                             const shared::linegraph::Route* r,
                                             const Node* dir,
                                             const Edge* fromEdge) const;

  bool containsRoute(const shared::linegraph::Route* r) const;
  size_t getCardinality() const;

  double getWidth() const;
  double getSpacing() const;

  double getTotalWidth() const;

  std::vector<const shared::linegraph::Route*> getShrdRoutes(
      const Edge& e) const;

  std::string toString() const;

 private:
  Node* _from;
  Node* _to;

  std::vector<RouteOccurance> _routes;

  PolyLine<double> _geom;

  double _width, _spacing;
};
}
}

#endif  // TRANSITMAP_GRAPH_EDGE_H_
