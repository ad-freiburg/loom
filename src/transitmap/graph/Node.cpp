// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <cassert>
#include "./Node.h"
#include "./Edge.h"
#include "./TransitGraph.h"
#include "gtfsparser/gtfs/Stop.h"
#include "../graph/EdgeTripGeom.h"

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
  if (s) _stops.insert(s);
}

// _____________________________________________________________________________
Node::Node(double x, double y, gtfs::Stop* s) : _pos(x, y) {
  if (s) _stops.insert(s);
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

// _____________________________________________________________________________
void Node::setPos(const util::geo::Point& p) {
  _pos = p;
}

// _____________________________________________________________________________
void Node::addMainDir(NodeFront f) {
  _mainDirs.push_back(f);
}

// _____________________________________________________________________________
const NodeFront* Node::getNodeFrontFor(const Edge* e) const {
  for (auto& nf : getMainDirs()) {
    if (std::find(nf.edges.begin(), nf.edges.end(), e) != nf.edges.end()) {
      return &nf;
    }
  }

  return 0;
}

// _____________________________________________________________________________
double Node::getScore() const {
  std::vector<geo::PolyLine> connections;
  for (auto nf : getMainDirs()) {
    for (auto e : nf.edges) {
      // TODO: only handle cases with 1 edgetripgeometry atm
      const EdgeTripGeom& etg = e->getEdgeTripGeoms()->front();
      size_t c = 0;
      for (auto rt : etg.getTrips()) {

      }
    }
  }

  return 0;
}

// _____________________________________________________________________________
std::vector<Partner> Node::getPartner(const NodeFront* f, const gtfs::Route* r) const {
  std::vector<Partner> ret;
  for (const auto& nf : getMainDirs()) {
    if (&nf == f) continue;

    for (const auto e : nf.edges) {
      for (const auto& etg : *e->getEdgeTripGeoms()) {
        for (const auto& to: etg.getTrips()) {
          if (to.route == r) {
            Partner p;
            p.front = &nf;
            p.edge = e;
            p.etg = &etg;
            p.route = to.route;
            ret.push_back(p);
          }
        }
      }
    }
  }
  return ret;
}
