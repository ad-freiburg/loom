// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/transitgraph/TransitEdgePL.h"
#include "shared/transitgraph/TransitGraph.h"
#include "shared/transitgraph/TransitNodePL.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdgePL;
using shared::transitgraph::RouteOcc;

// _____________________________________________________________________________
TransitEdgePL::TransitEdgePL() {}

// _____________________________________________________________________________
TransitEdgePL::TransitEdgePL(const PolyLine<double>& p)
    : _p(p) {}

// _____________________________________________________________________________
const util::geo::Line<double>* TransitEdgePL::getGeom() const {
  return &_p.getLine();
}

// _____________________________________________________________________________
const PolyLine<double>& TransitEdgePL::getPolyline() const { return _p; }

// _____________________________________________________________________________
void TransitEdgePL::setPolyline(const PolyLine<double>& p) { _p = p; }

// _____________________________________________________________________________
void TransitEdgePL::addRoute(
    const Route* r, const TransitNode* dir,
    util::Nullable<transitmapper::style::LineStyle> ls) {
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
void TransitEdgePL::addRoute(const Route* r, const TransitNode* dir) {
  addRoute(r, dir, util::Nullable<transitmapper::style::LineStyle>());
}

// _____________________________________________________________________________
void TransitEdgePL::delRoute(const Route* r) {
  RouteOcc occ(r, 0);
  _routes.erase(occ);
}

// _____________________________________________________________________________
const std::set<RouteOcc>& TransitEdgePL::getRoutes() const { return _routes; }

// _____________________________________________________________________________
std::set<RouteOcc>& TransitEdgePL::getRoutes() { return _routes; }

// _____________________________________________________________________________
util::json::Dict TransitEdgePL::getAttrs() const {
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
bool TransitEdgePL::hasRoute(const Route* r) const {
  return _routes.count(RouteOcc(r, 0)) > 0;
}

// _____________________________________________________________________________
const RouteOcc& TransitEdgePL::getRouteOcc(const Route* r) const {
  return *_routes.find(RouteOcc(r, 0));
}
