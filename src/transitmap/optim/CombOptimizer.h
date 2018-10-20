// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_COMBOPTIMIZER_H_
#define TRANSITMAP_OPTIM_COMBOPTIMIZER_H_

#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/Optimizer.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/NullOptimizer.h"
#include "transitmap/optim/Scorer.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

class CombOptimizer : public Optimizer {
 public:
  CombOptimizer(const config::Config* cfg, const Scorer* scorer)
      : _cfg(cfg), _scorer(scorer), _ilpOpt(cfg, scorer) {};

  int optimize(TransitGraph* tg) const;
  int optimize(const std::set<OptNode*>& g, HierarchOrderingConfig* c) const;

 private:
  const config::Config* _cfg;
  const Scorer* _scorer;
  const ILPEdgeOrderOptimizer _ilpOpt;
  const NullOptimizer _nullOpt;

  static size_t maxCard(const std::set<OptNode*>& g);
  static double solutionSpaceSize(const std::set<OptNode*>& g);

};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_COMBOPTIMIZER_H_
