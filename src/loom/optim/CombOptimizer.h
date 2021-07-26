// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_COMBOPTIMIZER_H_
#define LOOM_OPTIM_COMBOPTIMIZER_H_

#include "loom/config/TransitMapConfig.h"
#include "loom/graph/OrderCfg.h"
#include "loom/optim/ExhaustiveOptimizer.h"
#include "loom/optim/HillClimbOptimizer.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/Optimizer.h"
#include "loom/optim/Scorer.h"
#include "loom/optim/SimulatedAnnealingOptimizer.h"

using std::exception;
using std::string;

namespace loom {
namespace optim {

class CombOptimizer : public Optimizer {
 public:
  CombOptimizer(const config::Config* cfg, const Scorer* scorer)
      : Optimizer(cfg, scorer),
        _ilpOpt(cfg, scorer),
        _nullOpt(cfg, scorer),
        _exhausOpt(cfg, scorer),
        _hillcOpt(cfg, scorer),
        _annealOpt(cfg, scorer){};

  int optimizeComp(OptGraph* og, const std::set<OptNode*>& g, HierarOrderCfg* c,
                   size_t depth) const;

 private:
  const ILPEdgeOrderOptimizer _ilpOpt;
  const NullOptimizer _nullOpt;
  const ExhaustiveOptimizer _exhausOpt;
  const HillClimbOptimizer _hillcOpt;
  const SimulatedAnnealingOptimizer _annealOpt;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_COMBOPTIMIZER_H_
