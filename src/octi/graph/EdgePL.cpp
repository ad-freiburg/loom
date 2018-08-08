// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "octi/graph/EdgePL.h"
#include "octi/graph/NodePL.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using util::geo::PolyLine;
using namespace octi::graph;

// _____________________________________________________________________________
EdgePL::EdgePL() : _generation(-1) {}

// _____________________________________________________________________________
EdgePL::EdgePL(const PolyLine<double>& p) : _p(p), _generation(-1) {}

// _____________________________________________________________________________
const util::geo::Line<double>* EdgePL::getGeom() const { return &_p.getLine(); }

// _____________________________________________________________________________
const PolyLine<double>& EdgePL::getPolyline() const { return _p; }

// _____________________________________________________________________________
void EdgePL::setPolyline(const PolyLine<double>& p) { _p = p; }

// _____________________________________________________________________________
void EdgePL::addRoute(const Route* r, const Node<NodePL, EdgePL>* dir,
                      const LineStyle& ls) {
  _routes.insert(RouteOcc(r, dir, ls));
}

// _____________________________________________________________________________
void EdgePL::addRoute(const Route* r, const Node<NodePL, EdgePL>* dir) {
  _routes.insert(RouteOcc(r, dir));
}

// _____________________________________________________________________________
const std::set<RouteOcc>& EdgePL::getRoutes() const { return _routes; }

// _____________________________________________________________________________
void EdgePL::setGeneration(int64_t g) { _generation = g; }

// _____________________________________________________________________________
util::json::Dict EdgePL::getAttrs() const {
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
int64_t EdgePL::getGeneration() const { return _generation; }
