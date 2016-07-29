// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include "./node.h"
#include "./edge.h"
#include "./transitgraph.h"
#include "gtfsparser/gtfs/stop.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;

// _____________________________________________________________________________
Node::Node(util::geo::Point pos) : _pos(pos) {
}

// _____________________________________________________________________________
Node::Node(double x, double y) : _pos(x, y) {
}

// _____________________________________________________________________________
Node::Node(util::geo::Point pos, gtfs::Stop* s) : _pos(pos) {
  _stops.insert(s);
}

// _____________________________________________________________________________
Node::Node(double x, double y, gtfs::Stop* s) : _pos(x, y) {
  _stops.insert(s);
}

// _____________________________________________________________________________
Node::~Node() {
  for (Edge* e : _adjListOut) {
    e->getFrom()->removeEdge(e);
    e->getTo()->removeEdge(e);
    _adjListIn.erase(e); // catch edge to itself case
    delete e;
  }

  for (Edge* e : _adjListIn) {
    e->getFrom()->removeEdge(e);
    e->getTo()->removeEdge(e);
    delete e;
  }
}

// _____________________________________________________________________________
void Node::addStop(gtfs::Stop* s) {
  _stops.insert(s);
}

// _____________________________________________________________________________
const std::set<gtfs::Stop*>& Node::getStops() const {
  return _stops;
}

// _____________________________________________________________________________
void Node::addEdge(Edge* e) {
  if (e->getFrom() == this) _adjListOut.insert(e);
  if (e->getTo() == this) _adjListIn.insert(e);
}

// _____________________________________________________________________________
void Node::removeEdge(Edge* e) {
  if (e->getFrom() == this) _adjListOut.erase(e);
  if (e->getTo() == this) _adjListIn.erase(e);
}

// _____________________________________________________________________________
const util::geo::Point& Node::getPos() const {
  return _pos;
}
