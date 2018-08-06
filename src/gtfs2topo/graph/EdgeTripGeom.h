// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2TOPO_GRAPH_EDGETRIPGEOM_H_
#define GTFS2TOPO_GRAPH_EDGETRIPGEOM_H_

#include <vector>
#include "ad/cppgtfs/gtfs/Trip.h"
#include "ad/cppgtfs/gtfs/Route.h"
#include "util/geo/PolyLine.h"
#include "gtfs2topo/graph/BuildGraph.h"

namespace gtfs2topo {
namespace graph {

using namespace ad::cppgtfs;
using namespace util::geo;

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
  EdgeTripGeom(PolyLine<double> pl, const Node* geomDir);

  void addTrip(gtfs::Trip* t, const Node* dirNode, PolyLine<double>& pl);
  void addTrip(gtfs::Trip* t, const Node* dirNode);

  const std::vector<TripOccurance>& getTripsUnordered() const;
  std::vector<TripOccurance>* getTripsUnordered();

  std::vector<TripOccurance>::iterator removeTripOccurance(
      std::vector<TripOccurance>::const_iterator pos);

  TripOccurance* getTripsForRoute(const gtfs::Route* r) const;

  const PolyLine<double>& getGeom() const;
  void setGeom(const PolyLine<double>& p);

  bool containsRoute(gtfs::Route* r) const;
  size_t getTripCardinality() const;
  size_t getCardinality() const;
  const Node* getGeomDir() const;

  void setGeomDir(const Node* newDir);

  bool routeEquivalent(const EdgeTripGeom& g) const;

 private:
  std::vector<TripOccurance> _trips;

  PolyLine<double> _geom;
  const Node* _geomDir; // the direction of the geometry, may be reversed
};

}}

#endif  // GTFS2TOPO_GRAPH_EDGETRIPGEOM_H_

