// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_
#define TRANSITMAP_GRAPH_OPTIM_OPTGRAPH_H_

#include <proj_api.h>
#include <string>
#include <set>

#include "./../util/Geo.h"
#include "./OrderingConfiguration.h"
#include "./Node.h"
#include "./Edge.h"

using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::EdgeTripGeometry;

namespace transitmapper {
namespace optim {

struct OptNode {
  const Node* node;

  std::set<OptEdge*> adjListIn;
  std::set<OptEdge*> adjListOut;

  void Node::addEdge(OptEdge* e) {
    if (e->from == this) adjListOut.insert(e);
    if (e->to == this) adjListIn.insert(e);
  }
}

struct OptEdge {
  std::vector<EdgeTripGeometry*> etgs;

  const Node* from;
  const Node* to;

  OptEdge(const OptNode* from, const OptNode* to) : from(from), to(to);
}

class OptGraph {
 public:
  explicit OptGraph(TransitGraph* toOptim);

 private:
  TransitGraph* _g;
  std::set<OptNode*> _nodes;

  void build();
};

}}

#endif  // TRANSITMAP_GRAPH_TRANSITGRAPH_H_
