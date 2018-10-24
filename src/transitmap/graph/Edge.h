// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_EDGE_H_
#define TRANSITMAP_GRAPH_EDGE_H_

#include <vector>
#include "Edge.h"
#include "Node.h"
#include "Route.h"
#include "util/geo/PolyLine.h"
#include "util/Nullable.h"
#include "transitmap/style/LineStyle.h"

using std::exception;
using std::string;
using namespace util::geo;

namespace transitmapper {
namespace graph {

using style::LineStyle;

// forward declaration of Node
class Node;

struct RouteOccurance {
  RouteOccurance(const Route* r, const Node* dir) : route(r), direction(dir) {}
  RouteOccurance(const Route* r, const Node* dir, const style::LineStyle& ls)
   : route(r), direction(dir), style(ls) {}
  const Route* route;
  const Node* direction;  // 0 if in both directions

  util::Nullable<style::LineStyle> style;

  bool operator==(const RouteOccurance& b) const {
    return b.route == route;// && b.direction == direction;
  }
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

  void addRoute(const Route* r, const Node* dir, const LineStyle& ls);
  void addRoute(const Route* r, const Node* dir);

  const PolyLine<double>& getGeom() const;
  void setGeom(const PolyLine<double>& p);

  const std::vector<RouteOccurance>& getRoutes() const;
  std::vector<RouteOccurance>* getRoutes();

  RouteOccWithPos getRouteWithPosUnder(const Route* r,
    const std::vector<size_t> ordering) const;

  RouteOccWithPos getRouteWithPos(const Route* r) const;

  RouteOccurance* getRoute(const Route* r) const;

  std::set<const Route*> getRoutesRelTo(const Route* ref) const;

  std::vector<RouteOccurance> getCtdRoutesIn(const Node* n,
    const Route* r, const Node* dir, const Edge* fromEdge) const;

  bool containsRoute(const Route* r) const;
  size_t getCardinality() const;
  size_t getCardinality(bool woRelatives) const;

  double getWidth() const;
  double getSpacing() const;

  double getTotalWidth() const;

  std::vector<const Route*> getShrdRoutes(const Edge& e) const;

  bool replaceRoute(const Route* r, const Route* n);
  bool removeRoute(const Route* r);

  std::string toString() const;

  // TODO: store this here atm, but find better place...
  std::vector<std::vector<size_t> > permutations;
 private:
  Node* _from;
  Node* _to;

  std::vector<RouteOccurance> _routes;

  PolyLine<double> _geom;

  double _width, _spacing;
};

}}

#endif  // TRANSITMAP_GRAPH_EDGE_H_

