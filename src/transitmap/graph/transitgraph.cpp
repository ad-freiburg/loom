// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <string>
#include <set>
#include "transitgraph.h"
#include "edge.h"

using namespace transitmapper;
using namespace graph;

// _____________________________________________________________________________
TransitGraph::~TransitGraph() {
  // an edge is _deleted_ if either the from or to node is deleted!
  for (auto n : _nodes) {
    delete n;
  }
}


// _____________________________________________________________________________
TransitGraph::TransitGraph(const std::string& name) : _name(name) {

}

// _____________________________________________________________________________
Node* TransitGraph::addNode(Node* n) {
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
