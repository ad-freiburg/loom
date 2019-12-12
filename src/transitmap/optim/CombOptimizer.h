// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_COMBOPTIMIZER_H_
#define TRANSITMAP_OPTIM_COMBOPTIMIZER_H_

#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/ExhaustiveOptimizer.h"
#include "transitmap/optim/HillClimbOptimizer.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/NullOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/Optimizer.h"
#include "transitmap/optim/Scorer.h"
#include "transitmap/optim/SimulatedAnnealingOptimizer.h"

using std::exception;
using std::string;

namespace transitmapper {
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

  int optimizeComp(const std::set<OptNode*>& g, HierarchOrderingConfig* c) const;

 private:
  const ILPEdgeOrderOptimizer _ilpOpt;
  const NullOptimizer _nullOpt;
  const ExhaustiveOptimizer _exhausOpt;
  const HillClimbOptimizer _hillcOpt;
  const SimulatedAnnealingOptimizer _annealOpt;

};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_COMBOPTIMIZER_H_
