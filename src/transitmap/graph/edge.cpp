// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <vector>
#include "edge.h"
#include "node.h"
#include "gtfsparser/gtfs/trip.h"
#include "edgetripgeom.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;

// _____________________________________________________________________________
Edge::Edge(Node* from, Node* to) : _from(from), _to(to) {

}

// _____________________________________________________________________________
Node* Edge::getFrom() const {
  return _from;
}

// _____________________________________________________________________________
Node* Edge::getTo() const {
  return _to;
}

// _____________________________________________________________________________
bool Edge::addTrip(gtfs::Trip* t, Node* toNode) {
  for (auto& e : _tripsContained) {
    if (e.containsRoute(t->getRoute())) {
      for (auto& tr : (*e.getTrips())[t->getRoute()].trips) {
        // shorcut: if a trip is contained here with the same shape id,
        // don't require recalc of polyline etc
        if (tr->getShape() == t->getShape()) {
          e.addTrip(t, toNode);
          return true;
       }
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
bool Edge::addTrip(gtfs::Trip* t, geo::PolyLine pl, Node* toNode) {
  bool inserted = false;
  for (auto& e : _tripsContained) {
    if (e.getGeom() == pl || e.getGeom() == pl.getReversed()) {
      e.addTrip(t, toNode);
      inserted = true;
      break;
    }
  }

  if (!inserted) {
    EdgeTripGeom etg(pl);
    etg.addTrip(t, toNode);
    _tripsContained.push_back(etg);
  }

  return true;
}

// _____________________________________________________________________________
const std::vector<EdgeTripGeom>& Edge::getEdgeTripGeoms() const {
  return _tripsContained;
}

// _____________________________________________________________________________
std::vector<EdgeTripGeom>* Edge::getEdgeTripGeoms() {
  return &_tripsContained;
}

// _____________________________________________________________________________
void Edge::simplify() {
  for (auto& e : _tripsContained) {
    e.removeOrphans();
  }

// calculate average cardinalty of geometries on this edge
  double avg = 0;
  for (auto& e : _tripsContained) {
    avg += e.getTripCardinality();
  }

  avg /= _tripsContained.size();

  for (auto it = _tripsContained.begin(); it < _tripsContained.end(); ++it) {
    std::cout << it->getTripCardinality() << " vs. " << avg*0.1 << std::endl;
    if (it->getTripCardinality() < avg*0.1) {
      it = _tripsContained.erase(it);
    }
  }
}
