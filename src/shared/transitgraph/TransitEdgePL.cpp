// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/transitgraph/TransitEdgePL.h"
#include "shared/transitgraph/TransitNodePL.h"
#include "shared/transitgraph/TransitGraph.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdgePL;
using shared::transitgraph::RouteOcc;

// _____________________________________________________________________________
TransitEdgePL::TransitEdgePL() : _generation(-1) {}

// _____________________________________________________________________________
TransitEdgePL::TransitEdgePL(const PolyLine<double>& p)
    : _p(p), _generation(-1) {}

// _____________________________________________________________________________
const util::geo::Line<double>* TransitEdgePL::getGeom() const {
  return &_p.getLine();
}

// _____________________________________________________________________________
const PolyLine<double>& TransitEdgePL::getPolyline() const { return _p; }

// _____________________________________________________________________________
void TransitEdgePL::setPolyline(const PolyLine<double>& p) { _p = p; }

// _____________________________________________________________________________
void TransitEdgePL::addRoute(const Route* r, const TransitNode* dir,
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
    if (prev.direction != dir)  {
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
const std::set<RouteOcc>& TransitEdgePL::getRoutes() const { return _routes; }

// _____________________________________________________________________________
std::set<RouteOcc>& TransitEdgePL::getRoutes() { return _routes; }

// _____________________________________________________________________________
void TransitEdgePL::setGeneration(int64_t g) { _generation = g; }

// _____________________________________________________________________________
util::json::Dict TransitEdgePL::getAttrs() const {
  util::json::Dict obj;
  auto arr = util::json::Array();
  obj["generation"] = (int)_generation;

  for (auto r : getRoutes()) {
    auto route = util::json::Dict();
    route["id"] = r.route->getId();
    route["label"] = r.route->getLabel();
    route["color"] = r.route->getColor();

    if (r.direction != 0) {
      route["direction"] = util::toString(r.direction);
    }

    arr.push_back(route);
  }

  obj["lines"] = arr;

  return obj;
}

// _____________________________________________________________________________
bool TransitEdgePL::hasRoute(const Route* r) const {
  return _routes.count(RouteOcc(r, 0)) > 0;
}

// _____________________________________________________________________________
int64_t TransitEdgePL::getGeneration() const { return _generation; }
