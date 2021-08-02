// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_EXHAUSTIVEOPTIMIZER_H_
#define LOOM_OPTIM_EXHAUSTIVEOPTIMIZER_H_

#include "loom/config/LoomConfig.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "loom/optim/Optimizer.h"
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/RenderGraph.h"

namespace loom {
namespace optim {

class ExhaustiveOptimizer : public Optimizer {
 public:
  ExhaustiveOptimizer(const config::Config* cfg,
                      const shared::rendergraph::Penalties& pens)
      : Optimizer(cfg, pens), _optScorer(pens){};

  virtual int optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                           shared::rendergraph::HierarOrderCfg* c,
                           size_t depth) const;

 protected:
  OptGraphScorer _optScorer;
  void initialConfig(const std::set<OptNode*>& g, OptOrderCfg* cfg) const;
  void initialConfig(const std::set<OptNode*>& g, OptOrderCfg* cfg,
                     bool sorted) const;
  void writeHierarch(OptOrderCfg* cfg,
                     shared::rendergraph::HierarOrderCfg* c) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_EXHAUSTIVEOPTIMIZER_H_
