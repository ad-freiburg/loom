// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_EDGETRIPGEOM_H_
#define TRANSITMAP_GRAPH_EDGETRIPGEOM_H_

#include <vector>
#include "gtfsparser/gtfs/trip.h"
#include "gtfsparser/gtfs/route.h"
#include "../geo/PolyLine.h"

namespace transitmapper {
namespace graph {

using namespace gtfsparser;

class EdgeTripGeom {
 public:
  EdgeTripGeom(geo::PolyLine pl);
  void addTrip(gtfs::Trip* t);
  const std::map<gtfs::Route*, std::vector<gtfs::Trip*> >& getTrips() const;

  const geo::PolyLine& getGeom() const;
 private:
  std::map<gtfs::Route*, std::vector<gtfs::Trip*> > _trips;

  geo::PolyLine _geom;
};

}}

#endif  // TRANSITMAP_GRAPH_EDGE_H_

