// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include "ad/cppgtfs/gtfs/Stop.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "gtfs2graph/graph/EdgePL.h"
#include "gtfs2graph/graph/EdgeTripGeom.h"
#include "gtfs2graph/graph/NodePL.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/json/Writer.h"

using namespace gtfs2graph;
using namespace graph;
using namespace ad::cppgtfs;

using util::geo::DLine;
using util::geo::DPoint;

// _____________________________________________________________________________
NodePL::NodePL(DPoint pos) : _pos(pos) {}

// _____________________________________________________________________________
NodePL::NodePL(double x, double y) : _pos(x, y) {}

// _____________________________________________________________________________
NodePL::NodePL(DPoint pos, const gtfs::Stop* s) : _pos(pos) {
  if (s) _stops.insert(s);
}

// _____________________________________________________________________________
NodePL::NodePL(double x, double y, const gtfs::Stop* s) : _pos(x, y) {
  if (s) _stops.insert(s);
}

// _____________________________________________________________________________
void NodePL::addStop(const gtfs::Stop* s) { _stops.insert(s); }

// _____________________________________________________________________________
const std::set<const gtfs::Stop*>& NodePL::getStops() const { return _stops; }

// _____________________________________________________________________________
const DPoint& NodePL::getPos() const { return _pos; }

// _____________________________________________________________________________
void NodePL::setPos(const DPoint& p) { _pos = p; }

// _____________________________________________________________________________
bool NodePL::isConnOccuring(const gtfs::Route* r, const Edge* from,
                            const Edge* to) const {
  auto it = _occConns.find(r);
  if (it == _occConns.end()) return false;

  for (auto occ : it->second) {
    if ((occ.from == from && occ.to == to) ||
        (occ.from == to && occ.to == from)) {
      return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
void NodePL::connOccurs(const gtfs::Route* r, const Edge* from,
                        const Edge* to) {
  if (isConnOccuring(r, from, to)) return;

  _occConns[r].push_back(OccuringConnection(from, to));
}

// _____________________________________________________________________________
const std::map<const gtfs::Route*, std::vector<OccuringConnection> >&
NodePL::getOccuringConnections() const {
  return _occConns;
}

// _____________________________________________________________________________
std::vector<const Edge*> NodePL::getConnectingEdgesFor(const gtfs::Route* to,
                                                       Edge* a) const {
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
const util::geo::DPoint* NodePL::getGeom() const { return &getPos(); }

// _____________________________________________________________________________
void NodePL::setNode(const Node* n) { _n = n; }

// _____________________________________________________________________________
util::json::Dict NodePL::getAttrs() const {
  util::json::Dict obj;
  if (getStops().size() > 0) {
    obj["station_id"] = (*getStops().begin())->getId();
    obj["station_label"] = (*getStops().begin())->getName();
  }

  auto arr = util::json::Array();

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
            auto obj = util::json::Dict();
            obj["line"] = util::toString(r.route);
            obj["node_from"] =
                util::toString(e->getFrom() == _n ? e->getTo() : e->getFrom());
            obj["node_to"] =
                util::toString(f->getFrom() == _n ? f->getTo() : f->getFrom());
            arr.push_back(obj);
          }
        }
      }
    }
  }

  if (arr.size()) obj["excluded_conn"] = arr;
  return obj;
}
