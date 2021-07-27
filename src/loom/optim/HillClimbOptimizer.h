// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_HILLCLIMBOPTIMIZER_H_
#define LOOM_OPTIM_HILLCLIMBOPTIMIZER_H_

#include "loom/config/TransitMapConfig.h"
#include "loom/optim/ExhaustiveOptimizer.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "loom/optim/Optimizer.h"
#include "shared/rendergraph/OrderCfg.h"

namespace loom {
namespace optim {

class HillClimbOptimizer : public ExhaustiveOptimizer {
 public:
  HillClimbOptimizer(const config::Config* cfg, const Scorer* scorer)
      : ExhaustiveOptimizer(cfg, scorer){};

  virtual int optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                           shared::rendergraph::HierarOrderCfg* c,
                           size_t depth) const;

 private:
  double getScore(OptGraph* og, OptEdge* e, OptOrderCfg& cur) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_HILLCLIMBOPTIMIZER_H_
