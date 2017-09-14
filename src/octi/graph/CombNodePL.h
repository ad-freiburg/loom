// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_GRAPH_COMBNODEPL_H_
#define OCTI_GRAPH_COMBNODEPL_H_

#include "octi/graph/CombEdgePL.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"
#include "util/graph/Node.h"

using util::geo::Point;
using util::graph::Node;

namespace octi {
namespace graph {

typedef util::graph::Node<NodePL, EdgePL> Node;

// forward declaration
class CombNodePL;

struct PairCmp {
  bool operator()(const std::pair<util::graph::Edge<CombNodePL, CombEdgePL>*, double>& a,
                  const std::pair<util::graph::Edge<CombNodePL, CombEdgePL>*, double>& b) {
    return a.second < b.second;
  }
};

class CombNodePL : util::geograph::GeoNodePL {
 public:
  CombNodePL(octi::graph::Node* parent);

  const Point* getGeom() const;
  void getAttrs(json::object_t& obj) const;
  Node* getParent() const;

  void addOrderedEdge(util::graph::Edge<CombNodePL, CombEdgePL>* e, double deg);

  int32_t distBetween(util::graph::Edge<CombNodePL, CombEdgePL>* a, util::graph::Edge<CombNodePL, CombEdgePL>* b) const;

 private:
  octi::graph::Node* _parent;
  std::set<std::pair<util::graph::Edge<CombNodePL, CombEdgePL>*, double>, PairCmp>
      _edgeOrder;
};
}
}

#endif  // OCTI_GRAPH_COMBNODEPL_H_
