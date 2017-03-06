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

// _____________________________________________________________________________
const std::set<Edge*>& Node::getAdjListOut() const {
  return _adjListOut;
}

// _____________________________________________________________________________
const std::set<Edge*>& Node::getAdjListIn() const {
  return _adjListIn;
}

// _____________________________________________________________________________
std::set<Edge*> Node::getAdjList() const {
  std::set<Edge*> edges;
  edges.insert(getAdjListIn().begin(), getAdjListIn().end());
  edges.insert(getAdjListOut().begin(), getAdjListOut().end());

  return edges;
}

// _____________________________________________________________________________
bool Node::isConnOccuring(const gtfs::Route* r, const Edge* from, const Edge* to)
const {
  auto it = _occConns.find(r);
  if (it == _occConns.end()) return false;

  for (auto occ : it->second) {
    if ((occ.from == from && occ.to == to) || (occ.from == to && occ.to == from))
    {
      return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
void Node::connOccurs(const gtfs::Route* r, const Edge* from, const Edge* to) {
  if (isConnOccuring(r, from, to)) return;

  _occConns[r].push_back(OccuringConnection(from, to));
}

// _____________________________________________________________________________
void Node::replaceEdgeInConnections(const Edge* oldE, const Edge* newE) {
  for (auto it = _occConns.begin(); it != _occConns.end(); it++) {
    for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
      if (itt->from == oldE) itt->from = newE;
      if (itt->to == oldE) itt->to = newE;
    }
  }
}

// _____________________________________________________________________________
void Node::sewConnectionsTogether(const Edge* a, const Edge* b) {
  // TODO assertion that these edges are in here
  for (const auto& ega : a->getEdgeTripGeoms()) {
    for (const auto& to : ega.getTripsUnordered()) {
      for (const auto& egb : b->getEdgeTripGeoms()) {
        if (egb.containsRoute(to.route)) {
          std::cerr << "sweing..." << std::endl;
          connOccurs(to.route, a, b);
        }
      }
    }
  }
}