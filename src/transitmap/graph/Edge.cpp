// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "Edge.h"
#include "Node.h"
#include "util/geo/PolyLine.h"

using namespace transitmapper;
using namespace graph;

using util::geo::PolyLine;

// _____________________________________________________________________________
Edge::Edge(Node* from, Node* to, PolyLine pl, double w, double s)
    : _from(from), _to(to), _width(w), _spacing(s) {
  setGeom(pl);
}

// _____________________________________________________________________________
Node* Edge::getFrom() const { return _from; }

// _____________________________________________________________________________
Node* Edge::getTo() const { return _to; }

// _____________________________________________________________________________
void Edge::setFrom(Node* f) { _from = f; }

// _____________________________________________________________________________
void Edge::setTo(Node* t) { _to = t; }

// _____________________________________________________________________________
const std::vector<RouteOccurance>& Edge::getTripsUnordered() const {
  return _routes;
}

// _____________________________________________________________________________
Node* Edge::getOther(const Node* n) const {
  if (n == _from) return _to;
  if (n == _to) return _from;
  return 0;
}

// _____________________________________________________________________________
std::vector<RouteOccurance>* Edge::getTripsUnordered() { return &_routes; }

// _____________________________________________________________________________
std::vector<RouteOccurance> Edge::getContinuedRoutesIn(
    const Node* n, const Route* r, const Node* dir,
    const Edge* fromEdge) const {
  std::vector<RouteOccurance> ret;

  for (const RouteOccurance& to : _routes) {
    if (to.route == r) {
      if (to.direction == 0 || dir == 0 || (to.direction == n && dir != n) ||
          (to.direction != n && dir == n)) {
        if (n->connOccurs(r, fromEdge, this)) {
          ret.push_back(to);
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::set<const Route*> Edge::getRoutesRelTo(const Route* ref) const {
  std::set<const Route*> ret;

  for (const RouteOccurance& to : _routes) {
    if (to.route->relativeTo() == ref) {
      ret.insert(to.route);
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::vector<RouteOccurance> Edge::getSameDirRoutesIn(
    const Node* n, const Route* r, const Node* dir,
    const Edge* fromEdge) const {
  std::vector<RouteOccurance> ret;

  for (const RouteOccurance& to : _routes) {
    if (to.route == r) {

      if ((to.direction == 0 && dir == 0) ||
          (to.direction == n && dir != 0 && dir != n) ||
          (to.direction != n && to.direction != 0 && dir == n)) {
        if (n->connOccurs(r, fromEdge, this)) {
          ret.push_back(to);
        }
      } else {
        if (n->getId() == "0x2e138380") {
          std::cout << "fail" << std::endl;
          std::cout << to.route->getId() << " vs " << r->getId() << std::endl;
          std::cout << to.direction->getId() << " vs " << dir->getId() << std::endl;
          exit(0);
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
RouteOccurance* Edge::getRouteOcc(const Route* r) const {
  for (size_t i = 0; i < _routes.size(); i++) {
    RouteOccurance* to = const_cast<RouteOccurance*>(&_routes[i]);
    if (to->route == r) {
      return to;
    }
  }
  return 0;
}

// _____________________________________________________________________________
RouteOccWithPos Edge::getRouteOccWithPosUnder(
    const Route* r, const std::vector<size_t> ordering) const {
  for (size_t i = 0; i < _routes.size(); i++) {
    const RouteOccurance& to = _routes[i];
    if (to.route == r) {
      size_t pos =
          std::find(ordering.begin(), ordering.end(), i) - ordering.begin();
      return std::pair<RouteOccurance*, size_t>(
          const_cast<RouteOccurance*>(&to), pos);
    }
  }
  return std::pair<RouteOccurance*, size_t>(0, 0);
}

// _____________________________________________________________________________
RouteOccWithPos Edge::getRouteOccWithPos(const Route* r) const {
  for (size_t i = 0; i < _routes.size(); i++) {
    const RouteOccurance& to = _routes[i];
    if (to.route == r) {
      return std::pair<RouteOccurance*, size_t>(
          const_cast<RouteOccurance*>(&to), i);
    }
  }
  return std::pair<RouteOccurance*, size_t>(0, 0);
}

// _____________________________________________________________________________
void Edge::addRoute(const Route* r, const Node* dir, const LineStyle& ls) {
  if (containsRoute(r)) return;
  _routes.push_back(RouteOccurance(r, dir, ls));
}

// _____________________________________________________________________________
void Edge::addRoute(const Route* r, const Node* dir) {
  if (containsRoute(r)) return;
  _routes.push_back(RouteOccurance(r, dir));
}

// _____________________________________________________________________________
bool Edge::containsRoute(const Route* r) const {
  return getRouteOcc(r) != 0;
}

// _____________________________________________________________________________
size_t Edge::getCardinality() const { return _routes.size(); }

// _____________________________________________________________________________
size_t Edge::getCardinality(bool woRelatives) const {
  if (!woRelatives) return getCardinality();

  size_t ret = 0;

  for (const auto& ro : _routes) {
    if (ro.route->relativeTo() == 0) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
double Edge::getWidth() const { return _width; }

// _____________________________________________________________________________
double Edge::getSpacing() const { return _spacing; }

// _____________________________________________________________________________
double Edge::getTotalWidth() const {
  return getWidth() * _routes.size() + getSpacing() * (_routes.size() - 1);
}

// _____________________________________________________________________________
std::vector<const Route*> Edge::getSharedRoutes(const Edge& e) const {
  std::vector<const Route*> ret;
  for (auto& to : _routes) {
    if (e.containsRoute(to.route)) ret.push_back(to.route);
  }

  return ret;
}

// _____________________________________________________________________________
const PolyLine& Edge::getGeom() const { return _geom; }

// _____________________________________________________________________________
void Edge::setGeom(const PolyLine& p) { _geom = p; }
