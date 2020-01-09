// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/style/LineStyle.h"
#include "shared/linegraph/LineEdgePL.h"
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/LineNodePL.h"
#include "shared/linegraph/Route.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdgePL;
using shared::linegraph::RouteOcc;

// _____________________________________________________________________________
LineEdgePL::LineEdgePL() {}

// _____________________________________________________________________________
LineEdgePL::LineEdgePL(const PolyLine<double>& p) : _p(p) {}

// _____________________________________________________________________________
const util::geo::Line<double>* LineEdgePL::getGeom() const {
  return &_p.getLine();
}

// _____________________________________________________________________________
void LineEdgePL::setGeom(const util::geo::Line<double>& l) {
  _p = l;
}

// _____________________________________________________________________________
const PolyLine<double>& LineEdgePL::getPolyline() const { return _p; }

// _____________________________________________________________________________
void LineEdgePL::setPolyline(const PolyLine<double>& p) { _p = p; }

// _____________________________________________________________________________
void LineEdgePL::addRoute(const Route* r, const LineNode* dir,
                          util::Nullable<shared::style::LineStyle> ls) {
  RouteOcc occ(r, dir, ls);
  auto f = _routes.find(occ);
  if (f != _routes.end()) {
    const auto& prev = *f;
    // the route is already present in both directions, ignore newly inserted
    if (prev.direction == 0) return;

    // the route is already present in the same direction, ignore newly inserted
    if (prev.direction == dir) return;

    // the route is already present in the other direction, make two-way
    if (prev.direction != dir) {
      occ.direction = 0;
      _routes.erase(f);
    }
  }
  _routes.insert(occ);
}

// _____________________________________________________________________________
void LineEdgePL::addRoute(const Route* r, const LineNode* dir) {
  addRoute(r, dir, util::Nullable<shared::style::LineStyle>());
}

// _____________________________________________________________________________
void LineEdgePL::delRoute(const Route* r) {
  RouteOcc occ(r, 0);
  _routes.erase(occ);
}

// _____________________________________________________________________________
const std::set<RouteOcc>& LineEdgePL::getRoutes() const { return _routes; }

// _____________________________________________________________________________
std::set<RouteOcc>& LineEdgePL::getRoutes() { return _routes; }

// _____________________________________________________________________________
util::json::Dict LineEdgePL::getAttrs() const {
  util::json::Dict obj;
  auto arr = util::json::Array();

  std::string dbg_lines = "";
  bool first = true;

  for (auto r : getRoutes()) {
    auto route = util::json::Dict();
    route["id"] = r.route->getId();
    route["label"] = r.route->getLabel();
    route["color"] = r.route->getColor();

    if (r.direction != 0) {
      route["direction"] = util::toString(r.direction);
      dbg_lines += (first ? "" : "$") + r.route->getLabel() + ">";
    } else {
      dbg_lines += (first ? "" : "$") + r.route->getLabel();
    }

    arr.push_back(route);
    first = false;
  }

  obj["lines"] = arr;
  obj["dbg_lines"] = dbg_lines;

  return obj;
}

// _____________________________________________________________________________
bool LineEdgePL::hasRoute(const Route* r) const {
  return _routes.count(RouteOcc(r, 0)) > 0;
}

// _____________________________________________________________________________
const RouteOcc& LineEdgePL::getRouteOcc(const Route* r) const {
  return *_routes.find(RouteOcc(r, 0));
}

// _____________________________________________________________________________
const RouteOcc& LineEdgePL::routeOccAtPos(size_t i) const {
  auto it = _routes.begin();
  for (size_t j = 0; j < i; j++) it++;
  return *it;
}

// _____________________________________________________________________________
size_t LineEdgePL::getRoutePosUnder(
    const Route* r, const std::vector<size_t> ordering) const {
  size_t i = 0;
  for (const RouteOcc& ro : _routes) {
    if (ro.route == r) {
      size_t pos =
          std::find(ordering.begin(), ordering.end(), i) - ordering.begin();
      return pos;
    }
    i++;
  }
  return -1;
}

// _____________________________________________________________________________
size_t LineEdgePL::getRoutePos(const Route* r) const {
  size_t i = 0;
  for (const RouteOcc& ro : _routes) {
    if (ro.route == r) return i;
    i++;
  }
  return -1;
}
