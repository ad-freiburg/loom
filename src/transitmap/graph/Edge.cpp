// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <cassert>
#include <vector>
#include "Edge.h"
#include "Node.h"
#include "gtfsparser/gtfs/Trip.h"
#include "EdgeTripGeom.h"

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
  assert(toNode == _from || toNode == _to);
  for (auto& e : _tripsContained) {
    if (e.containsRoute(t->getRoute())) {
      for (auto& tr : e.getTripsForRoute(t->getRoute()).first->trips) {
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
void Edge::addEdgeTripGeom(const EdgeTripGeom& e) {
  assert(e.getGeomDir() == _from || e.getGeomDir() ==  _to);
  _tripsContained.push_back(e);
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
    if (it->getTripCardinality() < avg*0.1) {
      it = _tripsContained.erase(it);
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

  std::vector<const geo::PolyLine*> lines;
  std::vector<geo::PolyLine*> reversed;
  reversed.reserve(lines.size()); // prevent reallocation
  const Node* referenceDir = _tripsContained.front().getGeomDir();

  for (auto& et : _tripsContained) {
    if (et.getGeomDir() != referenceDir) {
      reversed.push_back(new geo::PolyLine(et.getGeom().getReversed()));
      lines.push_back(reversed.back());
    } else {
      lines.push_back(&et.getGeom());
    }
  }

  geo::PolyLine pl = geo::PolyLine::average(lines);

  EdgeTripGeom combined(pl, referenceDir);

  for (auto& et : _tripsContained) {
    for (auto& r : et.getTripsUnordered()) {
      for (auto& t : r.trips) {
        combined.addTrip(t, r.direction);
      }
    }
  }

  for (auto p : reversed) {
    delete p;
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
      if (toCheckAgainst.getGeom().contains(et->getGeom())
          && !et->getGeom().contains(toCheckAgainst.getGeom())) {
        for (auto& r : et->getTripsUnordered()) {
          for (auto& t : r.trips) {
            toCheckAgainst.addTrip(t, r.direction);
          }
        }
        combined = true;
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

// _____________________________________________________________________________
void Edge::fixEdgeTripGeomDirs() {
  for (auto& e : *getEdgeTripGeoms()) {
    if (e.getGeomDir() != _to) {
      e.setGeomDir(_to);
      (const_cast<geo::PolyLine*>(&e.getGeom()))->reverse();
    }
  }
}
