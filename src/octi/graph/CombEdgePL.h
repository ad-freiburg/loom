// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_COMBEDGEPL_H_
#define OCTI_GRAPH_COMBEDGEPL_H_

#include <set>
#include "util/geo/GeoGraph.h"
#include "octi/graph/Graph.h"

using util::geo::PolyLine;

namespace octi {
namespace graph {

typedef util::graph::Edge<NodePL, EdgePL> Edge;

class CombEdgePL : util::geograph::GeoEdgePL {
 public:
  CombEdgePL(Edge* child);

  const util::geo::Line* getGeom() const;
  void getAttrs(json::object_t& obj) const;

  std::vector<Edge*>& getChilds();

  const PolyLine& getPolyLine() const;
  void setPolyLine(const PolyLine& p);
  void setGeneration(int64_t g);
  int64_t getGeneration() const;
 private:
  std::vector<Edge*> _childs;

  int64_t _generation;

  PolyLine _geom;
};

}}

#endif  // OCTI_GRAPH_COMBEDGEPL_H_

