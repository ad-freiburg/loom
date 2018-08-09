// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_COMBEDGEPL_H_
#define OCTI_COMBGRAPH_COMBEDGEPL_H_

#include <set>
#include "octi/transitgraph/TransitGraph.h"
#include "util/geo/GeoGraph.h"

using util::geo::PolyLine;

namespace octi {
namespace combgraph {

class CombEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  CombEdgePL(octi::transitgraph::TransitEdge* child);

  const util::geo::Line<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  std::vector<octi::transitgraph::TransitEdge*>& getChilds();

  const PolyLine<double>& getPolyLine() const;
  void setPolyLine(const PolyLine<double>& p);
  void setGeneration(int64_t g);
  int64_t getGeneration() const;

 private:
  std::vector<octi::transitgraph::TransitEdge*> _childs;

  int64_t _generation;

  PolyLine<double> _geom;
};
}
}

#endif  // OCTI_COMBGRAPH_COMBEDGEPL_H_
