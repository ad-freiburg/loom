// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_EDGEPL_H_
#define OCTI_GRAPH_EDGEPL_H_

#include <set>
#include "util/geo/PolyLine.h"
#include "util/geo/GeoGraph.h"

using util::geo::PolyLine;

namespace octi {
namespace graph {

class EdgePL : util::geograph::GeoEdgePL {
 public:
  EdgePL(const PolyLine& p);

  void addRoute(const std::string& routeId);
  const std::set<std::string>& getRoutes() const;

  const util::geo::Line* getGeom() const;
  void getAttrs(json::object_t& obj) const;

  const PolyLine& getPolyline() const;
 private:
  std::set<std::string> _routes;

  PolyLine _p;
};

}}

#endif  // OCTI_GRAPH_EDGEPL_H_

