// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include <set>
#include <string>
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/OrderingConfiguration.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"

using util::geo::Point;
using util::geo::Box;
using util::geo::dist;
using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using transitmapper::graph::Route;
using transitmapper::graph::Configuration;

// _____________________________________________________________________________
TransitGraph::TransitGraph(const std::string& name, const std::string& proj)
    : _name(name) {
  _bbox = util::geo::minbox();
  _proj = pj_init_plus(proj.c_str());
}

// _____________________________________________________________________________
TransitGraph::~TransitGraph() {
  // an edge is _deleted_ if either the from or to node is deleted!
  // thus, we don't have to delete edges separately here
  for (auto n : _nodes) {
    delete n;
  }

  // clean up projection stuff
  pj_free(_proj);
}

// _____________________________________________________________________________
const std::string& TransitGraph::getName() const { return _name; }

// _____________________________________________________________________________
const Configuration& TransitGraph::getConfig() const { return _config; }

// _____________________________________________________________________________
void TransitGraph::setConfig(const Configuration& c) { _config = c; }

// _____________________________________________________________________________
double TransitGraph::getScore(double inStatPen, double sameSegCrossPen,
                              double diffSegCrossPen, double splitPen) const {
  return getScore(inStatPen, sameSegCrossPen, diffSegCrossPen, splitPen,
                  _config);
}

// _____________________________________________________________________________
double TransitGraph::getScore(double inStatPen, double sameSegCrossPen,
                              double diffSegCrossPen, double splitPen,
                              const Configuration& c) const {
  double ret = 0;

  for (auto n : getNodes()) {
    ret += n->getScore(inStatPen, sameSegCrossPen, diffSegCrossPen, splitPen,
                       true, true, c);
  }

  return ret;
}
//
// _____________________________________________________________________________
double TransitGraph::getCrossScore(double inStatPen, double sameSegCrossPen,
                                   double diffSegCrossPen) const {
  return getCrossScore(inStatPen, sameSegCrossPen, diffSegCrossPen, _config);
}

// _____________________________________________________________________________
double TransitGraph::getCrossScore(double inStatPen, double sameSegCrossPen,
                                   double diffSegCrossPen,
                                   const Configuration& c) const {
  double ret = 0;

  for (auto n : getNodes()) {
    ret += n->getCrossingScore(c, inStatPen, sameSegCrossPen, diffSegCrossPen,
                               true);
  }

  return ret;
}

// _____________________________________________________________________________
double TransitGraph::getSeparationScore(double inStatPen, double pen) const {
  return getSeparationScore(inStatPen, pen, _config);
}

// _____________________________________________________________________________
double TransitGraph::getSeparationScore(double inStatPen, double pen,
                                        const Configuration& c) const {
  double ret = 0;

  for (auto n : getNodes()) {
    ret += n->getSeparationScore(c, inStatPen, pen, true);
  }

  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumCrossings() const {
  return getNumCrossings(_config);
}

// _____________________________________________________________________________
size_t TransitGraph::getNumCrossings(const Configuration& c) const {
  double ret = 0;

  for (auto n : getNodes()) {
    ret += n->getNumCrossings(c);
  }

  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumSeparations() const {
  return getNumSeparations(_config);
}

// _____________________________________________________________________________
size_t TransitGraph::getNumSeparations(const Configuration& c) const {
  double ret = 0;

  for (auto n : getNodes()) {
    ret += n->getNumSeparations(c);
  }

  return ret;
}

// _____________________________________________________________________________
void TransitGraph::addNode(Node* n) {
  _nodes.insert(n);
  // expand the bounding box to hold this new node
  expandBBox(n->getPos());
}

// _____________________________________________________________________________
void TransitGraph::expandBBox(const Point& p) {
  bgeo::expand(_bbox, boost::geometry::make<bgeo::model::box<Point>>(
                          p.get<0>(), p.get<1>(), p.get<0>(), p.get<1>()));
}

// _____________________________________________________________________________
Node* TransitGraph::getNodeById(const std::string& id) const {
  for (auto n : _nodes) {
    if (n->getId() == id) return n;
  }

  return 0;
}

// _____________________________________________________________________________
Edge* TransitGraph::addEdge(Node* from, Node* to, PolyLine pl, double w,
                            double s) {
  if (from == to) return 0;
  Edge* e = getEdge(from, to);
  if (!e) {
    e = new Edge(from, to, pl, w, s);
    from->addEdge(e);
    to->addEdge(e);
    bgeo::expand(_bbox,
                 bgeo::return_envelope<bgeo::model::box<Point>>(pl.getLine()));
  }
  return e;
}

// _____________________________________________________________________________
void TransitGraph::deleteEdge(Node* from, Node* to) {
  Edge* toDel = getEdge(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdge(from, to));

  delete toDel;
}

// _____________________________________________________________________________
Route* TransitGraph::addRoute(const Route* r) {
  if (!getRoute(r->getId())) {
    _routes[r->getId()] = r;
  }
}

// _____________________________________________________________________________
const Route* TransitGraph::getRoute(const std::string& id) const {
  auto f = _routes.find(id);
  if (f == _routes.end()) return 0;

  return f->second;
}

// _____________________________________________________________________________
Edge* TransitGraph::getEdge(Node* from, Node* to) {
  for (auto e : from->getAdjListOut()) {
    if (e->getTo() == to) return e;
  }

  // also search in the opposite direction, we are handling an undirected
  // graph here
  for (auto e : from->getAdjListIn()) {
    if (e->getFrom() == to) return e;
  }

  return 0;
}

// _____________________________________________________________________________
const std::set<Node*>& TransitGraph::getNodes() const { return _nodes; }

// _____________________________________________________________________________
std::set<Node*>* TransitGraph::getNodes() { return &_nodes; }

// _____________________________________________________________________________
projPJ TransitGraph::getProjection() const { return _proj; }

// _____________________________________________________________________________
const Box& TransitGraph::getBoundingBox() const {
  return _bbox;
}

// _____________________________________________________________________________
Box TransitGraph::getBoundingBox(double p) const {
  return Box(
      Point(_bbox.min_corner().get<0>() - p, _bbox.min_corner().get<1>() - p),
      Point(_bbox.max_corner().get<0>() + p, _bbox.max_corner().get<1>() + p));
}

// _____________________________________________________________________________
Node* TransitGraph::getNearestNode(const Point& p, double maxD) const {
  double curD = DBL_MAX;
  ;
  Node* curN = 0;
  for (auto n : _nodes) {
    double d = dist(n->getPos(), p);
    if (d < maxD && d < curD) {
      curN = n;
      curD = d;
    }
  }

  return curN;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumNodes() const {
  return getNumNodes(true) + getNumNodes(false);
}

// _____________________________________________________________________________
size_t TransitGraph::getNumRoutes() const { return _routes.size(); }

// _____________________________________________________________________________
size_t TransitGraph::getMaxCardinality() const {
  size_t ret = 0;
  for (auto n : getNodes()) {
    for (auto e : n->getAdjListOut()) {
      if (e->getCardinality() > ret) ret = e->getCardinality();
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumEdges() const {
  size_t ret = 0;

  for (auto n : getNodes()) {
    ret += n->getAdjListOut().size();
  }

  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumNodes(bool topo) const {
  size_t ret = 0;
  for (auto n : _nodes) {
    if (n->getAdjListIn().size() + n->getAdjListOut().size() == 0) continue;
    if ((n->getStops().size() == 0) ^ !topo) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
double TransitGraph::getNumPossSolutions() const {
  double ret = 1;

  for (auto n : getNodes()) {
    for (auto e : n->getAdjListOut()) {
      ret *= util::factorial(e->getCardinality());
    }
  }

  return ret;
}

// _____________________________________________________________________________
double TransitGraph::getLastSolveTime() const {
  return _lastSolveTime;
}

// _____________________________________________________________________________
void TransitGraph::setLastSolveTime(double t) const {
  _lastSolveTime = t;
}
