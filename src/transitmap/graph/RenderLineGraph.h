// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_RENDERLINEGRAPH_H_
#define TRANSITMAP_GRAPH_RENDERLINEGRAPH_H_

#include <unordered_map>
#include "shared/transitgraph/TransitGraph.h"
#include "transitmap/config/TransitMapConfig.h"

using namespace util;

namespace transitmapper {
namespace graph {

using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdge;

class RenderLineGraph : public shared::transitgraph::TransitGraph {
 public:
  struct Front {
    Front(TransitEdge* edg, const util::geo::DLine& geom)
        : edg(edg), geom(geom){};
    TransitEdge* edg;
    util::geo::DLine geom;
  };

  RenderLineGraph(const transitmapper::config::Config* cfg) : _cfg(cfg){};
  ~RenderLineGraph(){};

  void addFront(const TransitNode* nd, const Front& fr);

  void writeStatGeoms();

 private:
  const transitmapper::config::Config* _cfg;
  std::unordered_map<const TransitNode*, std::vector<Front>> _fronts;
  std::unordered_map<const TransitNode*, util::geo::DPolygon> _hulls;

  Polygon<double> stationGeom(const TransitNode* node, double d,
                              bool rectangulize, bool simple) const;
};
}
}

#endif  // TRANSITMAP_GRAPH_RENDERLINEGRAPH_H_
