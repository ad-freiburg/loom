// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_EDGETRIPGEOM_H_
#define TRANSITMAP_GRAPH_EDGETRIPGEOM_H_

#include <vector>
#include "gtfsparser/gtfs/Trip.h"
#include "gtfsparser/gtfs/Route.h"
#include "../geo/PolyLine.h"
#include "./Node.h"

namespace transitmapper {
namespace graph {

using namespace gtfsparser;

struct TripOccurance {
  TripOccurance(gtfs::Route* r) : route(r), direction(0) {}
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
  gtfs::Route* route;
  std::vector<gtfs::Trip*> trips;
  const Node* direction;  // 0 if in both directions (should be 0 in most cases!)
};

typedef std::pair<TripOccurance*, size_t> TripOccWithPos;

class EdgeTripGeom {
 public:
  EdgeTripGeom(geo::PolyLine pl, const Node* geomDir);
  void addTrip(gtfs::Trip* t, const Node* dirNode, geo::PolyLine& pl);
  void addTrip(gtfs::Trip* t, const Node* dirNode);
  const std::vector<TripOccurance>& getTrips() const;
  std::vector<TripOccurance>* getTrips();

  TripOccWithPos getTripsForRoute(const gtfs::Route* r) const;

  const geo::PolyLine& getGeom() const;
  void setGeom(const geo::PolyLine& p);

  void removeOrphans();

  bool containsRoute(gtfs::Route* r) const;
  size_t getTripCardinality() const;
  const Node* getGeomDir() const;

  void setGeomDir(const Node* newDir);

  void setWidth(double w);
  void setSpacing(double s);

  double getWidth() const;
  double getSpacing() const;

  double getTotalWidth() const;
 private:
  std::vector<TripOccurance> _trips;

  double _w;
  double _s;
  geo::PolyLine _geom;
  const Node* _geomDir; // the direction of the geometry, may be reversed
};

}}

#endif  // TRANSITMAP_GRAPH_EDGE_H_

