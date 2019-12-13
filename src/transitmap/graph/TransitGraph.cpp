// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.  // Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <set>
#include <string>
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/OrderingConfig.h"
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
using transitmapper::graph::OrderingConfig;

// _____________________________________________________________________________
TransitGraph::TransitGraph()
     {
  _bbox = util::geo::minbox<double>();
}

// _____________________________________________________________________________
TransitGraph::~TransitGraph() {
  for (auto n : _nodes) {
    delete n;
  }
}

// _____________________________________________________________________________
const OrderingConfig& TransitGraph::getConfig() const { return _config; }

// _____________________________________________________________________________
void TransitGraph::setConfig(const OrderingConfig& c) { _config = c; }

// _____________________________________________________________________________
void TransitGraph::addNd(Node* n) {
  _nodes.insert(n);
  // expand the bounding box to hold this new node
  expandBBox(n->getPos());
}

// _____________________________________________________________________________
void TransitGraph::expandBBox(const DPoint& p) {
  _bbox = util::geo::extendBox(p, _bbox);
}

// _____________________________________________________________________________
Node* TransitGraph::getNodeById(const std::string& id) const {
  for (auto n : _nodes) {
    if (n->getId() == id) return n;
  }

  return 0;
}

// _____________________________________________________________________________
Edge* TransitGraph::addEdg(Node* from, Node* to, PolyLine<double> pl, double w,
                            double s) {
  if (from == to) return 0;
  Edge* e = getEdg(from, to);
  if (!e) {
    e = new Edge(from, to, pl, w, s);
    from->addEdg(e);
    to->addEdg(e);
    _bbox = util::geo::extendBox(util::geo::getBoundingBox(pl.getLine()), _bbox);
  }
  return e;
}

// _____________________________________________________________________________
void TransitGraph::deleteEdge(Node* from, Node* to) {
  Edge* toDel = getEdg(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdg(from, to));

  delete toDel;
}

// _____________________________________________________________________________
void TransitGraph::addRoute(const Route* r) {
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
Edge* TransitGraph::getEdg(Node* from, Node* to) {
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
const std::set<Node*>& TransitGraph::getNds() const { return _nodes; }

// _____________________________________________________________________________
std::set<Node*>* TransitGraph::getNds() { return &_nodes; }

// _____________________________________________________________________________
const DBox& TransitGraph::getBoundingBox() const {
  return _bbox;
}

// _____________________________________________________________________________
DBox TransitGraph::getBoundingBox(double p) const {
  return DBox(
      DPoint(_bbox.getLowerLeft().getX() - p, _bbox.getLowerLeft().getY() - p),
      DPoint(_bbox.getUpperRight().getX() + p, _bbox.getUpperRight().getY() + p));
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
  for (auto n : getNds()) {
    for (auto e : n->getAdjListOut()) {
      if (e->getCardinality() > ret) ret = e->getCardinality();
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getMaxDegree() const {
  size_t ret = 0;
  for (auto n : getNds()) {
    if (n->getAdjListOut().size() + n->getAdjListIn().size() > ret) {
      ret = n->getAdjListOut().size() + n->getAdjListIn().size();
    }
  }
  return ret;
}

// _____________________________________________________________________________
size_t TransitGraph::getNumEdges() const {
  size_t ret = 0;

  for (auto n : getNds()) {
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
