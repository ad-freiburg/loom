// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_NULLOPTIMIZER_H_
#define LOOM_OPTIM_NULLOPTIMIZER_H_

#include "loom/config/TransitMapConfig.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/Optimizer.h"
#include "loom/optim/Scorer.h"
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/RenderGraph.h"

namespace loom {
namespace optim {

class NullOptimizer : public Optimizer {
 public:
  NullOptimizer(const config::Config* cfg, const Scorer* scorer)
      : Optimizer(cfg, scorer){};
  int optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                   shared::rendergraph::HierarOrderCfg* c, size_t depth) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_NULLOPTIMIZER_H_
