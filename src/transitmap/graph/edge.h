// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_EDGE_H_
#define TRANSITMAP_GRAPH_EDGE_H_

#include <vector>
#include "edge.h"
#include "node.h"
#include "edgetripgeom.h"
#include "gtfsparser/gtfs/trip.h"

using namespace gtfsparser;

using std::exception;
using std::string;

namespace transitmapper {
namespace graph {

// forward declaration of Node
class Node;

class Edge {
 public:
  Edge(Node* from, Node* to);

  Node* getFrom() const;
  Node* getTo() const;

  void addTrip(gtfs::Trip* t);
 private:
  Node* _from;
  Node* _to;

  // Map of EdgeTripGeometries in this graph edge.
  // An EdgeTripGeometry is a geometry holding N trips.
  // This is meant as a multi-stage structure, where in the first
  // (trivial) stage, each EdgeTripGeometry holds exactly 1 trip.
  //
  // In a 2nd step, the EdgeTripGeometries are combined based on
  // geometrical equivalence.
  //
  // In a 3rd step, we split those edges with |_tripsContained| > 1
  // into single edges. This created either multiple distinct edges from
  // _from to _to, OR, in case the geometries have partial equivalence,
  // introduces new topological nodes which mark the position where two
  // lines part or join.
  std::vector<EdgeTripGeom> _tripsContained;

};

}}

#endif  // TRANSITMAP_GRAPH_EDGE_H_

