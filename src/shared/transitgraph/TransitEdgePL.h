// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_TRANSITGRAPH_TRANSITEDGEPL_H_
#define SHARED_TRANSITGRAPH_TRANSITEDGEPL_H_

#include <set>
#include "shared/style/LineStyle.h"
#include "shared/transitgraph/Route.h"
#include "util/Nullable.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Node.h"

namespace shared {
namespace transitgraph {

using util::geo::PolyLine;
using util::graph::Node;
using util::Nullable;

class TransitEdgePL;
class TransitNodePL;

struct RouteOcc {
  RouteOcc(const Route* r, const Node<TransitNodePL, TransitEdgePL>* dir)
      : route(r), direction(dir) {}
  RouteOcc(const Route* r, const Node<TransitNodePL, TransitEdgePL>* dir,
           const util::Nullable<shared::style::LineStyle>& ls)
      : route(r), direction(dir), style(ls) {}
  const Route* route;
  const Node<TransitNodePL, TransitEdgePL>*
      direction;  // 0 if in both directions

  util::Nullable<shared::style::LineStyle> style;
};

inline bool operator<(const RouteOcc& x, const RouteOcc& y) {
  return x.route < y.route;
};

class TransitEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  TransitEdgePL();
  TransitEdgePL(const PolyLine<double>& p);

  void addRoute(const Route* r, const Node<TransitNodePL, TransitEdgePL>* dir,
                util::Nullable<shared::style::LineStyle> ls);
  void addRoute(const Route* r, const Node<TransitNodePL, TransitEdgePL>* dir);

  const std::set<RouteOcc>& getRoutes() const;
  std::set<RouteOcc>& getRoutes();

  bool hasRoute(const Route* r) const;
  void delRoute(const Route* r);

  const RouteOcc& getRouteOcc(const Route* r) const;
  const RouteOcc& routeOccAtPos(size_t i) const;

  const util::geo::Line<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  const PolyLine<double>& getPolyline() const;
  void setPolyline(const PolyLine<double>& p);

 private:
  std::set<RouteOcc> _routes;

  PolyLine<double> _p;
};
}
}

#endif  // SHARED_TRANSITGRAPH_TRANSITEDGEPL_H_
