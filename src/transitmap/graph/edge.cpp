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
void Edge::addTrip(gtfs::Trip* t, geo::PolyLine pl) {
  bool inserted = false;
  for (auto& e : _tripsContained) {
    if (e.getGeom() == pl) {  // TODO: determine what equality means here
      e.addTrip(t);
      inserted = true;
      break;
    }
  }

  if (!inserted) {
    EdgeTripGeom etg(pl);
    etg.addTrip(t);
    _tripsContained.push_back(etg);
  }
}

// _____________________________________________________________________________
const std::vector<EdgeTripGeom>& Edge::getEdgeTripGeoms() const {
  return _tripsContained;
}
