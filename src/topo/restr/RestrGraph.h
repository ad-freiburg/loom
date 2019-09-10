// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_RESTR_RESTRGRAPH_H_
#define TOPO_RESTR_RESTRGRAPH_H_

#include <map>
#include <set>
#include "shared/transitgraph/TransitGraph.h"
#include "transitmap/graph/Route.h"
#include "util/graph/DirGraph.h"

using transitmapper::graph::Route;

namespace topo {
namespace restr {


struct RestrEdgePL;
struct RestrNodePL;

typedef util::graph::Node<RestrNodePL, RestrEdgePL> RestrNode;
typedef util::graph::Edge<RestrNodePL, RestrEdgePL> RestrEdge;

typedef std::map<const Route*,
                 std::map<const RestrEdge*, std::set<const RestrEdge*>>>
    ConnEx;

struct RestrNodePL {
  RestrNodePL() {};
  ConnEx restrs;
};

struct RestrEdgePL {
  RestrEdgePL(const util::geo::PolyLine<double>& geom) : geom(geom) {};
  util::geo::PolyLine<double> geom;
  std::set<Route*> routes;
};

typedef util::graph::DirGraph<RestrNodePL, RestrEdgePL> RestrGraph;
}
}

#endif  // TOPO_RESTR_RESTRINFERRER
