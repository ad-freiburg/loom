// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_RESTR_RESTRGRAPH_H_
#define TOPO_RESTR_RESTRGRAPH_H_

#include <map>
#include <set>
#include "shared/transitgraph/TransitGraph.h"
#include "shared/transitgraph/Route.h"
#include "util/graph/DirGraph.h"

namespace topo {
namespace restr {

struct RestrEdgePL;
struct RestrNodePL;

typedef util::graph::Node<RestrNodePL, RestrEdgePL> RestrNode;
typedef util::graph::Edge<RestrNodePL, RestrEdgePL> RestrEdge;

typedef std::map<const shared::transitgraph::Route*,
                 std::map<const RestrEdge*, std::set<const RestrEdge*>>>
    ConnEx;

struct RestrNodePL {
  RestrNodePL(){};
  ConnEx restrs;

  const util::geo::Point<double>* getGeom() const { return 0; }
  util::json::Dict getAttrs() const { return util::json::Dict(); }
};

struct RestrEdgePL {
  RestrEdgePL(const util::geo::PolyLine<double>& geom) : geom(geom){};
  util::geo::PolyLine<double> geom;
  std::set<const shared::transitgraph::Route*> routes;

  const util::geo::Line<double>* getGeom() const { return &geom.getLine(); }
  util::json::Dict getAttrs() const {
    util::json::Dict obj;
    auto arr = util::json::Array();

    std::string dbg_lines = "";
    bool first = true;

    for (auto r : routes) {
      auto route = util::json::Dict();
      route["id"] = r->getId();
      route["label"] = r->getLabel();
      route["color"] = r->getColor();

      dbg_lines += (first ? "" : "$") + r->getLabel();

      arr.push_back(route);
      first = false;
    }

    obj["lines"] = arr;
    obj["dbg_lines"] = dbg_lines;

    return obj;
  }
};

typedef util::graph::DirGraph<RestrNodePL, RestrEdgePL> RestrGraph;
}
}

#endif  // TOPO_RESTR_RESTRINFERRER
