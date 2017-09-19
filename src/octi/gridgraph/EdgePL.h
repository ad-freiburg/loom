// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRIDGRAPH_EDGEPL_H_
#define OCTI_GRIDGRAPH_EDGEPL_H_

#include <set>
#include "util/geo/PolyLine.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"
#include "octi/graph/CombEdgePL.h"
#include "octi/graph/CombNodePL.h"

using util::geo::PolyLine;

namespace octi {
namespace gridgraph {


class EdgePL : util::geograph::GeoEdgePL {
 public:
  EdgePL(const PolyLine& p, double c, bool secondar);

  void addRoute(std::string);
  const std::set<util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>*>& getResEdges() const;

  const util::geo::Line* getGeom() const;
  void getAttrs(json::object_t& obj) const;

  void setCost(double c);
  double cost() const;
  void addResidentEdge(util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>* e);
  bool isSecondary() const;
 private:
  PolyLine _pl;
  double _c;

  bool _isSecondary;

  std::set<util::graph::Edge<octi::graph::CombNodePL, octi::graph::CombEdgePL>*> _resEdges;
};

}}

#endif  // OCTI_GRIDGRAPH_EDGEPL_H_

