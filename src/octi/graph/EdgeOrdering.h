// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_EDGEORDERING_H_
#define OCTI_GRAPH_EDGEORDERING_H_

#include "octi/graph/CombEdgePL.h"
#include "util/graph/Node.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/Geo.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace graph {

// forward declaration
class CombNodePL;

typedef util::graph::Node<NodePL, EdgePL> Node;
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
  int64_t dist(CombEdge* a, CombEdge * b) const;
  const std::set<std::pair<CombEdge*, double>, PairCmp>& getOrderedSet() const;
  bool equals(const EdgeOrdering& e) const;
  std::string toString(CombNode* from) const;

 private:
  std::set<std::pair<CombEdge*, double>, PairCmp>
      _edgeOrder;
};

}}

#endif  // OCTI_GRAPH_EDGEORDERING_H_
