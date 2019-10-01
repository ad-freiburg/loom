// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_COMBGRAPH_EDGEORDERING_H_
#define OCTI_COMBGRAPH_EDGEORDERING_H_

#include <vector>
#include "octi/combgraph/CombEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;
using octi::combgraph::CombEdgePL;

namespace octi {
namespace combgraph {
// forward declaration
class CombNodePL;
}
}

using octi::combgraph::CombNodePL;

namespace octi {
namespace combgraph {

typedef util::graph::Node<CombNodePL, CombEdgePL> CombNode;
typedef util::graph::Edge<CombNodePL, CombEdgePL> CombEdge;

struct PairCmp {
  bool operator()(const std::pair<CombEdge*, double>& a,
                  const std::pair<CombEdge*, double>& b) {
    return a.second > b.second;
  }
};

class EdgeOrdering {
 public:
  void add(CombEdge* e, double deg);
  bool has(CombEdge* e) const;
  int64_t dist(CombEdge* a, CombEdge* b) const;
  const std::vector<std::pair<CombEdge*, double>>& getOrderedSet() const;
  bool equals(const EdgeOrdering& e) const;
  std::string toString(CombNode* from) const;

 private:
  std::vector<std::pair<CombEdge*, double>> _edgeOrder;
};
}
}

#endif  // OCTI_COMBGRAPH_EDGEORDERING_H_
