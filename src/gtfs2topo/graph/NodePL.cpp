// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include "gtfs2topo/graph/BuildGraph.h"
#include "util/geo/Geo.h"
#include "ad/cppgtfs/gtfs/Stop.h"
#include "gtfs2topo/graph/EdgeTripGeom.h"
#include "gtfs2topo/graph/NodePL.h"
#include "gtfs2topo/graph/EdgePL.h"
#include "util/String.h"

using namespace gtfs2topo;
using namespace graph;
using namespace ad::cppgtfs;

using util::geo::Point;
using util::geo::Line;


// _____________________________________________________________________________
NodePL::NodePL(Point pos) : _pos(pos) {
}

// _____________________________________________________________________________
NodePL::NodePL(double x, double y) : _pos(x, y) {
}

// _____________________________________________________________________________
NodePL::NodePL(Point pos, const gtfs::Stop* s) : _pos(pos) {
  if (s) _stops.insert(s);
}

// _____________________________________________________________________________
NodePL::NodePL(double x, double y, const gtfs::Stop* s) : _pos(x, y) {
  if (s) _stops.insert(s);
}

// _____________________________________________________________________________
void NodePL::addStop(const gtfs::Stop* s) {
  _stops.insert(s);
}

// _____________________________________________________________________________
const std::set<const gtfs::Stop*>& NodePL::getStops() const {
  return _stops;
}

// _____________________________________________________________________________
const Point& NodePL::getPos() const {
  return _pos;
}

// _____________________________________________________________________________
void NodePL::setPos(const Point& p) {
  _pos = p;
}

// _____________________________________________________________________________
bool NodePL::isConnOccuring(const gtfs::Route* r, const Edge* from, const Edge* to)
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
void NodePL::connOccurs(const gtfs::Route* r, const Edge* from, const Edge* to) {
  if (isConnOccuring(r, from, to)) return;

  _occConns[r].push_back(OccuringConnection(from, to));
}

// _____________________________________________________________________________
const std::map<const gtfs::Route*, std::vector<OccuringConnection> >& NodePL::getOccuringConnections()
const {
  return _occConns;
}

// _____________________________________________________________________________
void NodePL::replaceEdgeInConnections(const Edge* oldE, const Edge* newE) {
  for (auto it = _occConns.begin(); it != _occConns.end(); it++) {
    for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
      if (itt->from == oldE) itt->from = newE;
      if (itt->to == oldE) itt->to = newE;
    }
  }
}

// _____________________________________________________________________________
void NodePL::sewConnectionsTogether(Edge* a, Edge* b) {
  for (const auto& ega : *a->pl().getEdgeTripGeoms()) {
    for (const auto& to : ega.getTripsUnordered()) {
      for (const auto& egb : *b->pl().getEdgeTripGeoms()) {
        if (egb.containsRoute(to.route)) {
          connOccurs(to.route, a, b);
        }
      }
    }
  }
}

// _____________________________________________________________________________
std::vector<const Edge*> NodePL::getConnectingEdgesFor(const gtfs::Route* to, Edge* a)
const {
  std::vector<const Edge*> ret;

  auto it = _occConns.find(to);

  if (it == _occConns.end()) return ret;

  for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
    if (itt->from == a) ret.push_back(itt->to);
    if (itt->to == a) ret.push_back(itt->from);
  }

  return ret;
}

// _____________________________________________________________________________
const util::geo::Point* NodePL::getGeom() const {
  return &getPos();
}

// _____________________________________________________________________________
void NodePL::setNode(const Node* n) {
  _n = n;
}

// _____________________________________________________________________________
void NodePL::getAttrs(json::object_t& obj) const {
  if (getStops().size() > 0) {
    obj["station_id"] = (*getStops().begin())->getId();
    obj["station_label"] = (*getStops().begin())->getName();
  }

  auto arr = json::array();

  for (const graph::Edge* e : _n->getAdjList()) {
    if (!e->pl().getRefETG()) continue;
    for (auto r : e->pl().getRefETG()->getTripsUnordered()) {
      for (graph::Edge* f : _n->getAdjList()) {
        if (e == f) continue;
        if (!f->pl().getRefETG()) continue;
        for (auto rr : *f->pl().getRefETG()->getTripsUnordered()) {
          if (r.route == rr.route &&
              (r.direction == 0 || rr.direction == 0 ||
               (r.direction == _n && rr.direction != _n) ||
               (r.direction != _n && rr.direction == _n)) &&
              !isConnOccuring(r.route, e, f)) {
            auto obj = json::object();
            obj["route"] = util::toString(r.route);
            obj["edge1_node"] =
                util::toString(e->getFrom() == _n ? e->getTo() : e->getFrom());
            obj["edge2_node"] =
                util::toString(f->getFrom() == _n ? f->getTo() : f->getFrom());
            arr.push_back(obj);
          }
        }
      }
    }
  }

  if (arr.size()) obj["excluded_line_conns"] = arr;
}
