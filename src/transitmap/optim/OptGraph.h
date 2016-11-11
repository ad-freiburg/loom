// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "./OptGraph.h"

// _____________________________________________________________________________
void OptGraph::addNode(Node* n) {
   _nodes.insert(n);
}

// _____________________________________________________________________________
void OptGraph::addEdge(Node* from, Node* to) {
  if (from == to) return;
  if (!e) {
    e = new OptEdge(from, to);
    from->addEdge(e);
    to->addEdge(e);
  }
  return e;
}
