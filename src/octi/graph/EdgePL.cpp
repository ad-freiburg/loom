// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "util/geo/PolyLine.h"
#include "util/String.h"
#include "octi/graph/EdgePL.h"
#include "octi/graph/NodePL.h"

using util::geo::PolyLine;
using namespace octi::graph;

// _____________________________________________________________________________
EdgePL::EdgePL(const PolyLine& p) : _p(p), _generation(0) {

}

// _____________________________________________________________________________
const util::geo::Line* EdgePL::getGeom() const {
  return &_p.getLine();
}

// _____________________________________________________________________________
const PolyLine& EdgePL::getPolyline() const {
  return _p;
}

// _____________________________________________________________________________
void EdgePL::setPolyline(const PolyLine& p) {
  _p = p;
}

// _____________________________________________________________________________
void EdgePL::addRoute(const Route* r, const Node<NodePL, EdgePL>* dir, const LineStyle& ls) {
  _routes.insert(RouteOcc(r, dir, ls));
}

// _____________________________________________________________________________
void EdgePL::addRoute(const Route* r, const Node<NodePL, EdgePL>* dir) {
  _routes.insert(RouteOcc(r, dir));
}

// _____________________________________________________________________________
const std::set<RouteOcc>& EdgePL::getRoutes() const {
  return _routes;
}

// _____________________________________________________________________________
void EdgePL::setGeneration(size_t g) {
  _generation = g;
}

// _____________________________________________________________________________
void EdgePL::getAttrs(json::object_t& obj) const {
  obj["lines"] = json::array();
  obj["generation"] = _generation;

  for (auto r : getRoutes()) {
    json route = json::object();
    route["id"] = r.route->getId();
    route["label"] = r.route->getLabel();
    route["color"] = r.route->getColor();

    if (r.direction != 0) {
      route["direction"] = util::toString(r.direction);
    }

   obj["lines"].push_back(route);
  }
}
