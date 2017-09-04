// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_EDGEPL_H_
#define OCTI_GRIDGRAPH_EDGEPL_H_

#include "util/geo/PolyLine.h"
#include "util/geo/GeoGraph.h"

using util::geo::PolyLine;

namespace octi {
namespace gridgraph {

class EdgePL : util::geograph::GeoEdgePL {
 public:
  EdgePL(const PolyLine& p);

  const util::geo::Line* getGeom() const;
  void getAttrs(json::object_t& obj) const;
 private:
  PolyLine _pl;
};

}}

#endif  // OCTI_GRIDGRAPH_EDGEPL_H_

