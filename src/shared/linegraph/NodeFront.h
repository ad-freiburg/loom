// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_NODEFRONT_H_
#define SHARED_LINEGRAPH_NODEFRONT_H_

#include "shared/linegraph/LineGraph.h"
#include "shared/linegraph/LineEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"
#include "util/graph/Edge.h"
#include "transitmap/graph/OrderingConfig.h"

namespace shared {
namespace linegraph {

struct NodeFront {
  NodeFront(LineNode* n, LineEdge* e) : n(n), edge(e) {}

  util::geo::DPoint getTripOccPos(const shared::linegraph::Route* r,
                       const transitmapper::graph::OrderingConfig& c, bool origGeom) const;
  util::geo::DPoint getTripPos(const LineEdge* e, size_t pos, bool inv, bool originGeom) const;

  double getOutAngle() const;

  LineNode* n;
  LineEdge* edge;

  // geometry after expansion
  PolyLine<double> geom;

  // geometry before expansion
  PolyLine<double> origGeom;

  void setInitialGeom(const PolyLine<double>& g) {
    geom = g;
    origGeom = g;
  };
  void setGeom(const PolyLine<double>& g) { geom = g; };

  // TODO
  double refEtgLengthBefExp;
};
}
}

#endif  // SHARED_LINEGRAPH_NODEFRONT_H_
