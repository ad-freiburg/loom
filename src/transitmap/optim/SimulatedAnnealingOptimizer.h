// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_SIMULATEDANNEALINGOPTIMIZER_H_
#define TRANSITMAP_OPTIM_SIMULATEDANNEALINGOPTIMIZER_H_

#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/NullOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/OptGraphScorer.h"
#include "transitmap/optim/Optimizer.h"
#include "transitmap/optim/HillClimbOptimizer.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

class SimulatedAnnealingOptimizer : public HillClimbOptimizer {
 public:
  SimulatedAnnealingOptimizer(const config::Config* cfg, const Scorer* scorer)
      : HillClimbOptimizer(cfg, scorer) {};

  virtual int optimize(const std::set<OptNode*>& g, HierarchOrderingConfig* c) const;

 private:
  double getScore(OptEdge* e, OptOrderingConfig& cur) const;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_SIMULATEDANNEALINGOPTIMIZER_H_
