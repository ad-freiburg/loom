// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_EDGEPL_H_
#define OCTI_GRAPH_EDGEPL_H_

#include <set>
#include "util/geo/PolyLine.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"
#include "transitmap/graph/Route.h"
#include "transitmap/style/LineStyle.h"
#include "util/Nullable.h"

using util::geo::PolyLine;
using transitmapper::graph::Route;
using transitmapper::style::LineStyle;
using util::graph::Node;
using util::Nullable;

namespace octi {
namespace graph {

class EdgePL;
class NodePL;

struct RouteOcc {
  RouteOcc(const Route* r, const Node<NodePL, EdgePL>* dir) : route(r), direction(dir) {}
  RouteOcc(const Route* r, const Node<NodePL, EdgePL>* dir, const transitmapper::style::LineStyle& ls)
   : route(r), direction(dir), style(ls) {}
  const Route* route;
  const Node<NodePL, EdgePL>* direction;  // 0 if in both directions

  util::Nullable<transitmapper::style::LineStyle> style;
};

inline bool operator<(const RouteOcc& x, const RouteOcc& y) {
  return x.route < y.route;
};

class EdgePL : util::geograph::GeoEdgePL {
 public:
  EdgePL(const PolyLine& p);

  void addRoute(const Route* r, const Node<NodePL, EdgePL>* dir, const LineStyle& ls);
  void addRoute(const Route* r, const Node<NodePL, EdgePL>* dir);

  const std::set<RouteOcc>& getRoutes() const;

  const util::geo::Line* getGeom() const;
  void getAttrs(json::object_t& obj) const;

  const PolyLine& getPolyline() const;
  void setPolyline(const PolyLine& p);
 private:
  std::set<RouteOcc> _routes;

  PolyLine _p;
};

}}

#endif  // OCTI_GRAPH_EDGEPL_H_

