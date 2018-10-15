// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/OptGraph.h"

#ifndef TRANSITMAP_OPTIM_OPTIMIZER_H_
#define TRANSITMAP_OPTIM_OPTIMIZER_H_

using transitmapper::graph::OrderingConfig;
using transitmapper::graph::HierarchOrderingConfig;
using transitmapper::graph::Route;
using transitmapper::graph::TransitGraph;
using transitmapper::optim::OptNode;
using transitmapper::optim::OptEdge;

namespace transitmapper {
namespace optim {

class Optimizer {
 public:
  virtual int optimize(TransitGraph* tg) const = 0;
  virtual int optimize(const std::set<OptNode*>& g,
                       HierarchOrderingConfig* c) const = 0;

 protected:
  static void expandRelatives(TransitGraph* g, OrderingConfig* c);
  static void expandRelativesFor(
      OrderingConfig* c, const Route* ref, graph::Edge* start,
      const std::set<const Route*>& rs,
      std::set<std::pair<graph::Edge*, const Route*>>& visited);
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_OPTIMIZER_H_
