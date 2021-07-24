// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_EXHAUSTIVEOPTIMIZER_H_
#define LOOM_OPTIM_EXHAUSTIVEOPTIMIZER_H_

#include "loom/config/TransitMapConfig.h"
#include "loom/graph/OrderCfg.h"
#include "loom/graph/RenderGraph.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "loom/optim/Optimizer.h"

using std::exception;
using std::string;

namespace loom {
namespace optim {

class ExhaustiveOptimizer : public Optimizer {
 public:
  ExhaustiveOptimizer(const config::Config* cfg, const Scorer* scorer)
      : Optimizer(cfg, scorer), _optScorer(scorer){};

  virtual int optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                           HierarOrderCfg* c, size_t depth) const;

 protected:
  OptGraphScorer _optScorer;
  void initialConfig(const std::set<OptNode*>& g, OptOrderCfg* cfg) const;
  void initialConfig(const std::set<OptNode*>& g, OptOrderCfg* cfg,
                     bool sorted) const;
  void writeHierarch(OptOrderCfg* cfg, HierarOrderCfg* c) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_EXHAUSTIVEOPTIMIZER_H_
