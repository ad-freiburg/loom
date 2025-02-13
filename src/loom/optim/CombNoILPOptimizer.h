// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_COMBNOILPOPTIMIZER_H_
#define LOOM_OPTIM_COMBNOILPOPTIMIZER_H_

#include "loom/config/LoomConfig.h"
#include "loom/optim/ExhaustiveOptimizer.h"
#include "loom/optim/HillClimbOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/Optimizer.h"
#include "loom/optim/SimulatedAnnealingOptimizer.h"
#include "shared/rendergraph/OrderCfg.h"

namespace loom {
namespace optim {

class CombNoILPOptimizer : public Optimizer {
 public:
  CombNoILPOptimizer(const config::Config* cfg,
                     const shared::rendergraph::Penalties& pens)
      : Optimizer(cfg, pens),
        _nullOpt(cfg, pens),
        _exhausOpt(cfg, pens),
        _hillcOpt(cfg, pens, false),
        _annealOpt(cfg, pens, false){};

  double optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                      shared::rendergraph::HierarOrderCfg* c, size_t depth,
                      OptResStats& stats) const;

  virtual std::string getName() const { return "comb-no-ilp"; }

 private:
  const NullOptimizer _nullOpt;
  const ExhaustiveOptimizer _exhausOpt;
  const HillClimbOptimizer _hillcOpt;
  const SimulatedAnnealingOptimizer _annealOpt;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_COMBNOILPOPTIMIZER_H_
