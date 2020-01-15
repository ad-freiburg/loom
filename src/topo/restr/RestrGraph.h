// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_RESTR_RESTRGRAPH_H_
#define TOPO_RESTR_RESTRGRAPH_H_

#include <map>
#include <set>
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/Line.h"
#include "util/graph/DirGraph.h"

namespace topo {
namespace restr {

struct RestrEdgePL;
struct RestrNodePL;

typedef util::graph::Node<RestrNodePL, RestrEdgePL> RestrNode;
typedef util::graph::Edge<RestrNodePL, RestrEdgePL> RestrEdge;

typedef std::map<const shared::linegraph::Line*,
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
  std::set<const shared::linegraph::Line*> lines;

  const util::geo::Line<double>* getGeom() const { return &geom.getLine(); }
  util::json::Dict getAttrs() const {
    util::json::Dict obj;
    auto arr = util::json::Array();

    std::string dbg_lines = "";
    bool first = true;

    for (auto l : lines) {
      auto line = util::json::Dict();
      line["id"] = l->id();
      line["label"] = l->label();
      line["color"] = l->color();

      dbg_lines += (first ? "" : "$") + l->label();

      arr.push_back(line);
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
