// Copyright 2022, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPOEVAL_DIRLINEGRAPH_H_
#define TOPOEVAL_DIRLINEGRAPH_H_

#include <map>
#include <set>
#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/Line.h"
#include "util/graph/DirGraph.h"

namespace topoeval {

struct DirLineEdgePL;
struct DirLineNodePL;

typedef util::graph::Node<DirLineNodePL, DirLineEdgePL> DirLineNode;
typedef util::graph::Edge<DirLineNodePL, DirLineEdgePL> DirLineEdge;

typedef std::map<const shared::linegraph::Line*,
                 std::map<const DirLineEdge*, std::set<const DirLineEdge*>>>
    ConnEx;

struct DirLineNodePL {
  DirLineNodePL() {};
  DirLineNodePL(const std::string& stationLabel, const util::geo::DPoint& geom) : _statLbl(stationLabel), _geom(geom) {};
  ConnEx restrs;

  const std::string& getStatLabel() const { return _statLbl; }
  const util::geo::Point<double>* getGeom() const { return &_geom; }
  util::json::Dict getAttrs() const;

  std::string _statLbl;
  util::geo::DPoint _geom;
};

struct DirLineEdgePL {
  DirLineEdgePL(const util::geo::PolyLine<double>& geom) : geom(geom){};
  util::geo::PolyLine<double> geom;
  std::set<const shared::linegraph::Line*> lines;

  const util::geo::Line<double>* getGeom() const { return 0; }
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

typedef util::graph::DirGraph<DirLineNodePL, DirLineEdgePL> DirLineGraph;

// _____________________________________________________________________________
inline util::json::Dict DirLineNodePL::getAttrs() const {
  util::json::Dict obj;
  auto arr = util::json::Array();

  for (const auto& ro : restrs) {
    for (const auto& exFr : ro.second) {
      for (const auto* exTo : exFr.second) {
        util::json::Dict ex;
        ex["line"] = util::toString(ro.first->id());
        if (exFr.first == exTo) continue;
        auto shrd = DirLineGraph::sharedNode(exFr.first, exTo);
        auto nd1 = exFr.first->getOtherNd(shrd);
        auto nd2 = exTo->getOtherNd(shrd);
        ex["node_from"] = util::toString(nd1);
        ex["node_to"] = util::toString(nd2);
        arr.push_back(ex);
      }
    }
  }

  if (arr.size()) obj["excluded_conn"] = arr;
  return obj;
}
}

#endif  // TOPOEVAL_DIRLINEGRAPH
