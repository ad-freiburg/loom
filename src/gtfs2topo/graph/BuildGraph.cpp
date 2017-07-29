// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include <set>
#include <string>
#include "util/geo/Geo.h"
#include "gtfs2topo/graph/Edge.h"
#include "gtfs2topo/graph/BuildGraph.h"

using gtfs2topo::graph::BuildGraph;
using gtfs2topo::graph::Node;
using gtfs2topo::graph::Edge;
using util::geo::Point;

// _____________________________________________________________________________
BuildGraph::BuildGraph(const std::string& name, const std::string& proj) : _name(name) {
  _proj = pj_init_plus(proj.c_str());
}

// _____________________________________________________________________________
BuildGraph::~BuildGraph() {
  // an edge is _deleted_ if either the from or to node is deleted!
  // thus, we don't have to delete edges separately here
  for (auto n : _nodes) {
    delete n;
  }

  // clean up projection stuff
  pj_free(_proj);
}

// _____________________________________________________________________________
Node* BuildGraph::addNode(Node* n) {
  auto ins = _nodes.insert(n);
  return *ins.first;
}

// _____________________________________________________________________________
Edge* BuildGraph::addEdge(Node* from, Node* to) {
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
Edge* BuildGraph::getEdge(Node* from, Node* to) {
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
const std::set<Node*>& BuildGraph::getNodes() const { return _nodes; }

// _____________________________________________________________________________
std::set<Node*>* BuildGraph::getNodes() { return &_nodes; }

// _____________________________________________________________________________
Node* BuildGraph::getNodeByStop(const gtfs::Stop* s, bool getParent) const {
  if (getParent && s->getParentStation())
    return getNodeByStop(s->getParentStation());

  return getNodeByStop(s);
}

// _____________________________________________________________________________
Node* BuildGraph::getNodeByStop(const gtfs::Stop* s) const {
  for (const auto n : _nodes) {
    if (n->getStops().find(const_cast<gtfs::Stop*>(s)) != n->getStops().end()) {
      return n;
    }
  }
  return 0;
}

// _____________________________________________________________________________
void BuildGraph::deleteEdge(Node* from, Node* to) {
  Edge* toDel = getEdge(from, to);
  if (!toDel) return;

  from->removeEdge(toDel);
  to->removeEdge(toDel);

  assert(!getEdge(from, to));

  delete toDel;
}

// _____________________________________________________________________________
void BuildGraph::deleteNode(Node* n) { _nodes.erase(n); }

// _____________________________________________________________________________
projPJ BuildGraph::getProjection() const { return _proj; }

// _____________________________________________________________________________
Node* BuildGraph::getNearestNode(const Point& p, double maxD) const {
  double curD = DBL_MAX;
  ;
  Node* curN = 0;
  for (auto n : _nodes) {
    double d = util::geo::dist(n->getPos(), p);
    if (d < maxD && d < curD) {
      curN = n;
      curD = d;
    }
  }

  return curN;
}
