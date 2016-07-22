// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <proj_api.h>
#include <string>
#include <set>
#include "transitgraph.h"
#include "edge.h"

using namespace transitmapper;
using namespace graph;


// _____________________________________________________________________________
TransitGraph::TransitGraph(const std::string& name, const std::string& proj)
: _name(name) {
  _bbox = boost::geometry::make_inverse<boost::geometry::model::box<util::geo::Point> >();
  _proj = pj_init_plus(proj.c_str());
}

// _____________________________________________________________________________
TransitGraph::~TransitGraph() {
  // an edge is _deleted_ if either the from or to node is deleted!
  for (auto n : _nodes) {
    delete n;
  }

  // clean up projection stuff
  pj_free(_proj);
}

// _____________________________________________________________________________
Node* TransitGraph::addNode(Node* n) {
  // expand the bounding box to hold this new node
  boost::geometry::expand(_bbox, n->getPos());

  std::cout << boost::geometry::wkt(n->getPos()) << std::endl;
  std::cout << boost::geometry::wkt(_bbox) << std::endl;
  return *(_nodes.insert(n)).first;
}

// _____________________________________________________________________________
Edge* TransitGraph::addEdge(Node* from, Node* to) {
  Edge* e = new Edge(from, to);
  from->addEdge(e);
  to->addEdge(e);
  return e;
}

// _____________________________________________________________________________
const std::set<Node*>& TransitGraph::getNodes() const {
  return _nodes;
}

// _____________________________________________________________________________
std::set<Node*>* TransitGraph::getNodes() {
  return &_nodes;
}

// _____________________________________________________________________________
bool TransitGraph::containsNode(Node* n) const {
  return  _nodes.find(n) != _nodes.end();
}

// _____________________________________________________________________________
void TransitGraph::deleteNode(Node* n) {
  _nodes.erase(n);
}

// _____________________________________________________________________________
projPJ TransitGraph::getProjection() const {
  return _proj;
}

// _____________________________________________________________________________
const boost::geometry::model::box<util::geo::Point>& TransitGraph::getBoundingBox() const {
  return _bbox;
}
