// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_EXHAUSTIVEOPTIMIZER_H_
#define TRANSITMAP_OPTIM_EXHAUSTIVEOPTIMIZER_H_

#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/NullOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/OptGraphScorer.h"
#include "transitmap/optim/Optimizer.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

class ExhaustiveOptimizer : public Optimizer {
 public:
  ExhaustiveOptimizer(const config::Config* cfg, const Scorer* scorer)
      : Optimizer(cfg, scorer), _optScorer(scorer){};

  virtual int optimizeComp(const std::set<OptNode*>& g,
                           HierarchOrderingConfig* c, size_t depth) const;

 protected:
  OptGraphScorer _optScorer;
  void initialConfig(const std::set<OptNode*>& g, OptOrderingConfig* cfg) const;
  void initialConfig(const std::set<OptNode*>& g, OptOrderingConfig* cfg,
                     bool sorted) const;
  void writeHierarch(OptOrderingConfig* cfg, HierarchOrderingConfig* c) const;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_EXHAUSTIVEOPTIMIZER_H_
