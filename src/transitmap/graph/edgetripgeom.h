// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_EDGETRIPGEOM_H_
#define TRANSITMAP_GRAPH_EDGETRIPGEOM_H_

#include <vector>
#include "gtfsparser/gtfs/trip.h"
#include "gtfsparser/gtfs/route.h"
#include "../geo/PolyLine.h"
#include "./node.h"

namespace transitmapper {
namespace graph {

using namespace gtfsparser;

struct TripOccurance {
  TripOccurance() : direction(0) {}
  void addTrip(gtfs::Trip* t, const Node* dirNode) {
    if (trips.size() == 0) {
      direction = dirNode;
    } else {
      if (direction && direction != dirNode) {
        direction = 0;
      }
    }
    trips.push_back(t);
  }
  std::vector<gtfs::Trip*> trips;
  const Node* direction;  // 0 if in both directions (should be 0 in most cases!)
};

class EdgeTripGeom {
 public:
  EdgeTripGeom(geo::PolyLine pl, const Node* geomDir);
  void addTrip(gtfs::Trip* t, const Node* dirNode, geo::PolyLine& pl);
  void addTrip(gtfs::Trip* t, const Node* dirNode);
  const std::map<gtfs::Route*, TripOccurance>& getTrips() const;
  std::map<gtfs::Route*, TripOccurance>* getTrips();

  const geo::PolyLine& getGeom() const;
  void setGeom(const geo::PolyLine& p);

  void removeOrphans();

  bool containsRoute(gtfs::Route* r) const;
  size_t getTripCardinality() const;
  const Node* getGeomDir() const;

  void setGeomDir(const Node* newDir);
 private:
  std::map<gtfs::Route*, TripOccurance> _trips;

  geo::PolyLine _geom;
  const Node* _geomDir; // the direction of the geometry, may be reversed
};

}}

#endif  // TRANSITMAP_GRAPH_EDGE_H_

