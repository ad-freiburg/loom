// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include <string>
#include <set>
#include "./Graph.h"
#include "./Edge.h"

using skeletonbuilder::graph::Graph;
using skeletonbuilder::graph::Node;
using skeletonbuilder::graph::Edge;
using pbutil::geo::Point;
using bgeo::make_inverse;

// _____________________________________________________________________________
Graph::Graph(const std::string& name, const std::string& proj)
: _name(name) {
  _proj = pj_init_plus(proj.c_str());
}

// _____________________________________________________________________________
Graph::~Graph() {
  // an edge is _deleted_ if either the from or to node is deleted!
  // thus, we don't have to delete edges separately here
  for (auto n : _nodes) {
    delete n;
  }

  // clean up projection stuff
  pj_free(_proj);
}

// _____________________________________________________________________________
Node* Graph::addNode(Node* n) {
  auto ins = _nodes.insert(n);
  return *ins.first;
}

// _____________________________________________________________________________
Edge* Graph::addEdge(Node* from, Node* to) {
  if (from == to) return 0;
  Edge* e = getEdge(from, to);
  if (!e) {
    e = new Edge(from, to);
    from->addEdge(e);
    to->addEdge(e);
  }
  return e;
}

// _____________________________________________________________________________
Edge* Graph::getEdge(Node* from, Node* to) {
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
const std::set<Node*>& Graph::getNodes() const {
  return _nodes;
}

// _____________________________________________________________________________
std::set<Node*>* Graph::getNodes() {
  return &_nodes;
}

// _____________________________________________________________________________
Node* Graph::getNodeByStop(const gtfs::Stop* s, bool getParent) const {
  if (getParent && s->getParentStation()) return getNodeByStop(
    s->getParentStation());

  return getNodeByStop(s);
}

// _____________________________________________________________________________
Node* Graph::getNodeByStop(const gtfs::Stop* s) const {
  for (const auto n : _nodes) {
    if (n->getStops().find(const_cast<gtfs::Stop*>(s)) != n->getStops().end()) {
      return n;
    }
  }
  return 0;
}

// _____________________________________________________________________________
void Graph::deleteEdge(Node* from, Node* to) {
  Edge* toDel = getEdge(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdge(from, to));

  delete toDel;
}

// _____________________________________________________________________________
void Graph::deleteNode(Node* n) {
  _nodes.erase(n);
}

// _____________________________________________________________________________
projPJ Graph::getProjection() const {
  return _proj;
}

// _____________________________________________________________________________
Node* Graph::getNearestNode(const Point& p, double maxD) const {
  double curD = DBL_MAX;;
  Node* curN = 0;
  for (auto n : _nodes) {
    double d = bgeo::distance(n->getPos(), p);
    if (d < maxD && d < curD) {
      curN = n;
      curD = d;
    }
  }

  return curN;
}
