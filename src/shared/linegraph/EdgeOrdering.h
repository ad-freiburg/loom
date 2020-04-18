// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_EDGEORDERING_H_
#define SHARED_LINEGRAPH_EDGEORDERING_H_

#include <vector>
#include "shared/linegraph/LineEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;
using shared::linegraph::LineEdgePL;

namespace shared {
namespace linegraph {

typedef util::graph::Node<LineNodePL, LineEdgePL> LineNode;
typedef util::graph::Edge<LineNodePL, LineEdgePL> LineEdge;

struct PairCmp {
  bool operator()(const std::pair<LineEdge*, double>& a,
                  const std::pair<LineEdge*, double>& b) {
    return a.second > b.second;
  }
};

class EdgeOrdering {
 public:
  void add(LineEdge* e, double deg);
  bool has(LineEdge* e) const;
  int64_t dist(LineEdge* a, LineEdge* b) const;
  const std::vector<std::pair<LineEdge*, double>>& getOrderedSet() const;
  bool equals(const EdgeOrdering& e) const;
  std::string toString(LineEdge* from) const;

 private:
  std::vector<std::pair<LineEdge*, double>> _edgeOrder;
};
}
}

#endif  // SHARED_LINEGRAPH_EDGEORDERING_H_
