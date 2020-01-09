// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINEEDGEPL_H_
#define SHARED_LINEGRAPH_LINEEDGEPL_H_

#include <set>
#include "shared/style/LineStyle.h"
#include "shared/linegraph/Route.h"
#include "util/Nullable.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Node.h"

namespace shared {
namespace linegraph {

using util::geo::PolyLine;
using util::graph::Node;
using util::Nullable;

class LineEdgePL;
class LineNodePL;

struct RouteOcc {
  RouteOcc(const Route* r, const Node<LineNodePL, LineEdgePL>* dir)
      : route(r), direction(dir) {}
  RouteOcc(const Route* r, const Node<LineNodePL, LineEdgePL>* dir,
           const util::Nullable<shared::style::LineStyle>& ls)
      : route(r), direction(dir), style(ls) {}
  const Route* route;
  const Node<LineNodePL, LineEdgePL>*
      direction;  // 0 if in both directions

  util::Nullable<shared::style::LineStyle> style;
};

inline bool operator<(const RouteOcc& x, const RouteOcc& y) {
  return x.route < y.route;
};

class LineEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  LineEdgePL();
  LineEdgePL(const PolyLine<double>& p);

  void addRoute(const Route* r, const Node<LineNodePL, LineEdgePL>* dir,
                util::Nullable<shared::style::LineStyle> ls);
  void addRoute(const Route* r, const Node<LineNodePL, LineEdgePL>* dir);

  const std::set<RouteOcc>& getRoutes() const;
  std::set<RouteOcc>& getRoutes();

  bool hasRoute(const Route* r) const;
  void delRoute(const Route* r);

  const RouteOcc& getRouteOcc(const Route* r) const;
  const RouteOcc& routeOccAtPos(size_t i) const;

  size_t getRoutePosUnder(const Route* r,
                          const std::vector<size_t> ordering) const;

  size_t getRoutePos(const Route* r) const;

  const util::geo::Line<double>* getGeom() const;
  void setGeom(const util::geo::Line<double>& l);
  util::json::Dict getAttrs() const;

  const PolyLine<double>& getPolyline() const;
  void setPolyline(const PolyLine<double>& p);

 private:
  std::set<RouteOcc> _routes;

  PolyLine<double> _p;
};
}
}

#endif  // SHARED_LINEGRAPH_LINEEDGEPL_H_
