// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef GTFS2GRAPH_GRAPH_EDGEPL_H_
#define GTFS2GRAPH_GRAPH_EDGEPL_H_

#include <vector>
#include "ad/cppgtfs/gtfs/Trip.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "gtfs2graph/graph/EdgeTripGeom.h"
#include "util/Nullable.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/PolyLine.h"
#include "util/json/Writer.h"

using namespace ad::cppgtfs;
using util::Nullable;

namespace gtfs2graph {
namespace graph {

class EdgePL : public util::geograph::GeoEdgePL<double> {
 public:
  EdgePL(const Edge* e);
  EdgePL();

  bool addTrip(gtfs::Trip* t, util::geo::PolyLine<double> pl, Node* toNode);

  void setEdge(const Edge* e);

  const std::vector<EdgeTripGeom>& getEdgeTripGeoms() const;
  std::vector<EdgeTripGeom>* getEdgeTripGeoms();

  const EdgeTripGeom* getRefETG() const;
  EdgeTripGeom* getRefETG();

  void addEdgeTripGeom(const EdgeTripGeom& e);

  void simplify(double pruneThreshold);
  void setTo(Node* to);

  const util::geo::Line<double>* getGeom() const;
  util::json::Dict getAttrs() const;

 private:
  std::vector<EdgeTripGeom> _tripsContained;

  void combineIncludedGeoms();
  void averageCombineGeom();

  const Edge* _e;
};

}  // namespace graph
}  // namespace gtfs2graph

#endif  // GTFS2GRAPH_GRAPH_EDGEPL_H_
