// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_ORACLEOPTIMIZER_H_
#define LOOM_OPTIM_ORACLEOPTIMIZER_H_

#include "loom/config/LoomConfig.h"
#include "loom/optim/ExhaustiveOptimizer.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/NullOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "shared/rendergraph/OrderCfg.h"
#include "shared/rendergraph/RenderGraph.h"

namespace loom {
namespace optim {

class OracleOptimizer : public ExhaustiveOptimizer {
 public:
  OracleOptimizer(const config::Config* cfg,
                  const shared::rendergraph::Penalties& pens)
      : ExhaustiveOptimizer(cfg, pens){};

  virtual int optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                           shared::rendergraph::HierarOrderCfg* c,
                           size_t depth) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_ORACLEOPTIMIZER_H_
