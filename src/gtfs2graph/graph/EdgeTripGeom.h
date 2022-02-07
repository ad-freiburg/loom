// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2GRAPH_GRAPH_EDGETRIPGEOM_H_
#define GTFS2GRAPH_GRAPH_EDGETRIPGEOM_H_

#include <vector>
#include "ad/cppgtfs/gtfs/Route.h"
#include "ad/cppgtfs/gtfs/Trip.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "util/geo/PolyLine.h"

namespace gtfs2graph {
namespace graph {

struct RouteOccurance {
  RouteOccurance(ad::cppgtfs::gtfs::Route* r) : route(r), direction(0) {}
  void addTrip(ad::cppgtfs::gtfs::Trip* t, const Node* dirNode) {
    if (trips.size() == 0) {
      direction = dirNode;
    } else {
      if (direction && direction != dirNode) direction = 0;
    }
    trips.push_back(t);
  }
  ad::cppgtfs::gtfs::Route* route;
  std::vector<ad::cppgtfs::gtfs::Trip*> trips;
  const Node* direction;  // 0 if in both directions
};

typedef std::pair<RouteOccurance*, size_t> TripOccWithPos;

class EdgeTripGeom {
 public:
  EdgeTripGeom(util::geo::PolyLine<double> pl, const Node* geomDir);

  void addTrip(ad::cppgtfs::gtfs::Trip* t, const Node* dirNode,
               util::geo::PolyLine<double>& pl);
  void addTrip(ad::cppgtfs::gtfs::Trip* t, const Node* dirNode);

  const std::vector<RouteOccurance>& getTripsUnordered() const;
  std::vector<RouteOccurance>* getTripsUnordered();

  std::vector<RouteOccurance>::iterator removeRouteOcc(
      std::vector<RouteOccurance>::const_iterator pos);

  RouteOccurance* getRouteOcc(const ad::cppgtfs::gtfs::Route* r) const;

  const util::geo::PolyLine<double>& getGeom() const;
  void setGeom(const util::geo::PolyLine<double>& p);

  bool containsRoute(ad::cppgtfs::gtfs::Route* r) const;
  size_t getTripCardinality() const;
  size_t getCardinality() const;
  const Node* getGeomDir() const;

  void setGeomDir(const Node* newDir);

 private:
  std::vector<RouteOccurance> _routeOccs;

  util::geo::PolyLine<double> _geom;
  const Node* _geomDir;  // the direction of the geometry, may be reversed
};

}  // namespace graph
}  // namespace gtfs2graph

#endif  // GTFS2GRAPH_GRAPH_EDGETRIPGEOM_H_
