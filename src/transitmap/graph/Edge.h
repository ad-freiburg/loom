// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_EDGE_H_
#define TRANSITMAP_GRAPH_EDGE_H_

#include <vector>
#include "Edge.h"
#include "Node.h"
#include "Route.h"
#include "pbutil/geo/PolyLine.h"
#include "pbutil/Nullable.h"
#include "../style/LineStyle.h"

using std::exception;
using std::string;
using namespace pbutil::geo;

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

  pbutil::Nullable<style::LineStyle> style;
};

typedef std::pair<RouteOccurance*, size_t> RouteOccWithPos;

class Edge {
 public:
  Edge(Node* from, Node* to, PolyLine pl, double w, double s);

  Node* getFrom() const;
  Node* getTo() const;

  void setFrom(Node* from);
  void setTo(Node* to);

  void addRoute(const Route* r, const Node* dir, const LineStyle& ls);
  void addRoute(const Route* r, const Node* dir);

  const PolyLine& getGeom() const;
  void setGeom(const PolyLine& p);

  const std::vector<RouteOccurance>& getTripsUnordered() const;
  std::vector<RouteOccurance>* getTripsUnordered();

  RouteOccWithPos getTripsForRouteUnder(const Route* r,
    const std::vector<size_t> ordering) const;

  RouteOccurance* getTripsForRoute(const Route* r) const;

  std::vector<RouteOccurance> getContinuedRoutesIn(const Node* n,
    const Route* r, const Node* dir, const Edge* fromEdge) const;

  std::vector<RouteOccurance> getSameDirRoutesIn(const Node* n,
    const Route* r, const Node* dir, const Edge* fromEdge) const;

  bool containsRoute(const Route* r) const;
  size_t getCardinality() const;

  double getWidth() const;
  double getSpacing() const;

  double getTotalWidth() const;

  std::vector<const Route*> getSharedRoutes(const Edge& e) const;


  // TODO: store this here atm, but find better plcae...
  std::vector<std::vector<size_t> > permutations;
 private:
  Node* _from;
  Node* _to;

  std::vector<RouteOccurance> _routes;

  PolyLine _geom;

  double _width, _spacing;
};

}}

#endif  // TRANSITMAP_GRAPH_EDGE_H_

