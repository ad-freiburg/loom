// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_COMBEDGEPL_H_
#define OCTI_COMBGRAPH_COMBEDGEPL_H_

#include <set>
#include "shared/transitgraph/TransitGraph.h"
#include "util/geo/GeoGraph.h"

using util::geo::PolyLine;

namespace octi {
namespace combgraph {

class CombEdgePL : util::geograph::GeoEdgePL<double> {
 public:
  CombEdgePL(shared::transitgraph::TransitEdge* child);

  const util::geo::Line<double>* getGeom() const;
  util::json::Dict getAttrs() const;

  const std::vector<shared::transitgraph::TransitEdge*>& getChilds() const;
  std::vector<shared::transitgraph::TransitEdge*>& getChilds();

  const PolyLine<double>& getPolyLine() const;
  void setPolyLine(const PolyLine<double>& p);

 private:
  std::vector<shared::transitgraph::TransitEdge*> _childs;

  PolyLine<double> _geom;
};
}
}

#endif  // OCTI_COMBGRAPH_COMBEDGEPL_H_
