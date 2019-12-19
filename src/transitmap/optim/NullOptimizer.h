// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_NULLOPTIMIZER_H_
#define TRANSITMAP_OPTIM_NULLOPTIMIZER_H_

#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/Optimizer.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/NullOptimizer.h"
#include "transitmap/optim/Scorer.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

class NullOptimizer : public Optimizer {
 public:
  NullOptimizer(const config::Config* cfg, const Scorer* scorer)
      : Optimizer(cfg, scorer) {};
  int optimizeComp(const std::set<OptNode*>& g, HierarchOrderingConfig* c, size_t depth) const;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_NULLOPTIMIZER_H_
