// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "gtfsparser/gtfs/Trip.h"
#include "pbutil/geo/PolyLine.h"
#include "Edge.h"
#include "Node.h"
#include "EdgeTripGeom.h"

using namespace skeletonbuilder;
using namespace graph;
using namespace gtfsparser;

using pbutil::geo::PolyLine;

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
Node* Edge::getOtherNode(const Node* notNode) const {
  if (_to == notNode) return _from;
  return _to;
}

// _____________________________________________________________________________
bool Edge::addTrip(gtfs::Trip* t, Node* toNode) {
  assert(toNode == _from || toNode == _to);
  for (auto& e : _tripsContained) {
    if (e.containsRoute(t->getRoute())) {
      for (auto& tr : e.getTripsForRoute(t->getRoute())->trips) {
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
bool Edge::addTrip(gtfs::Trip* t, PolyLine pl, Node* toNode) {
  assert(toNode == _from || toNode == _to);
  bool inserted = false;
  for (auto& e : _tripsContained) {
    if (e.getGeom() == pl) {
      e.addTrip(t, toNode, pl);
      inserted = true;
      break;
    }
  }
  if (!inserted) {
    EdgeTripGeom etg(pl, toNode);
    etg.addTrip(t, toNode);
    addEdgeTripGeom(etg);
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
void Edge::addEdgeTripGeom(const EdgeTripGeom& e) {
  assert(e.getGeomDir() == _from || e.getGeomDir() ==  _to);

  _tripsContained.push_back(e);
  if (e.getGeomDir() != _to) {
    const_cast<PolyLine*>(&_tripsContained.back().getGeom())->reverse();
    _tripsContained.back().setGeomDir(_to);
  }
}

// _____________________________________________________________________________
void Edge::simplify() {
  /**
  for (auto& e : _tripsContained) {
    e.removeOrphans();
  }
  **/

  // calculate average cardinalty of geometries on this edge
  double avg = 0;
  for (auto& e : _tripsContained) {
    avg += e.getTripCardinality();
  }

  avg /= _tripsContained.size();

  for (auto it = _tripsContained.begin(); it < _tripsContained.end(); ++it) {
    if (it->getTripCardinality() < avg*0.1) {
      //it = _tripsContained.erase(it);
    }
  }

  combineIncludedGeoms();
  averageCombineGeom();
}

// _____________________________________________________________________________
void Edge::averageCombineGeom() {
  if (_tripsContained.size() < 2) {
    return;
  }

  std::vector<const PolyLine*> lines;

  for (auto& et : _tripsContained) {
    assert(et.getGeomDir() == _to);
    lines.push_back(&et.getGeom());
  }

  PolyLine pl = PolyLine::average(lines);

  EdgeTripGeom combined(pl, _to);

  for (auto& et : _tripsContained) {
    for (auto& r : *et.getTripsUnordered()) {
      for (auto& t : r.trips) {
        combined.addTrip(t, r.direction);
      }
    }
  }

  _tripsContained.clear();
  _tripsContained.push_back(combined);
}

// _____________________________________________________________________________
void Edge::combineIncludedGeoms() {
  if (_tripsContained.size() < 2) {
    return;
  }

  for (auto et = _tripsContained.begin(); et != _tripsContained.end();) {
    bool combined = false;
    for (auto& toCheckAgainst : _tripsContained) {
      if (toCheckAgainst.getGeom().getLength() > et->getGeom().getLength()
          && toCheckAgainst.getGeom().contains(et->getGeom(), 50)
          && !et->getGeom().contains(toCheckAgainst.getGeom(), 50)) {
        for (auto& r : *et->getTripsUnordered()) {
          for (auto& t : r.trips) {
            toCheckAgainst.addTrip(t, r.direction);
          }
        }
        combined = true;
        break;
      }
    }
    if (combined) {
      // delete the old EdgeTripGeom
      et = _tripsContained.erase(et);
    } else {
      et++;
    }
  }
}
