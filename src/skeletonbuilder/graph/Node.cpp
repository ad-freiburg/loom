// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include "./Node.h"
#include "./Edge.h"
#include "./Graph.h"
#include "transitmap/util/Geo.h"
#include "gtfsparser/gtfs/Stop.h"
#include "./EdgeTripGeom.h"

using namespace transitmapper;
using namespace skeletonbuilder;
using namespace graph;
using namespace gtfsparser;

using util::geo::Point;
using util::geo::Line;


// _____________________________________________________________________________
Node::Node(Point pos) : _pos(pos) {
}

// _____________________________________________________________________________
Node::Node(double x, double y) : _pos(x, y) {
}

// _____________________________________________________________________________
Node::Node(Point pos, gtfs::Stop* s) : _pos(pos) {
  if (s) _stops.insert(s);
}

// _____________________________________________________________________________
Node::Node(double x, double y, gtfs::Stop* s) : _pos(x, y) {
  if (s) _stops.insert(s);
}

// _____________________________________________________________________________
Node::~Node() {
  for (auto e = _adjListOut.begin(); e != _adjListOut.end();) {
    Edge* eP = *e;

    if (eP->getFrom() == this) {
      // careful with invalidating iterators
      e = _adjListOut.erase(e);
    } else {
      eP->getFrom()->removeEdge(eP);
      e++;
    }

    eP->getTo()->removeEdge(eP);

    delete eP;
  }

  for (auto e = _adjListIn.begin(); e != _adjListIn.end();) {
    Edge* eP = *e;

    if (eP->getTo() == this) {
      // careful with invalidating iterators
      e = _adjListIn.erase(e);
    } else {
      eP->getTo()->removeEdge(eP);
      e++;
    }

    eP->getFrom()->removeEdge(eP);

    delete eP;
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
const Point& Node::getPos() const {
  return _pos;
}

// _____________________________________________________________________________
void Node::setPos(const Point& p) {
  _pos = p;
}