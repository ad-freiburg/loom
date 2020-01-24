// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_
#define TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_

#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/OrderCfg.h"
#include "transitmap/optim/ILPOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/Optimizer.h"
#include "transitmap/optim/Scorer.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace optim {

using namespace graph;

typedef std::pair<size_t, size_t> PosCom;
typedef std::pair<PosCom, PosCom> PosComPair;
typedef std::pair<OptEdge*, OptEdge*> EdgePair;

class ILPEdgeOrderOptimizer : public ILPOptimizer {
 public:
  ILPEdgeOrderOptimizer(const config::Config* cfg, const Scorer* scorer)
      : ILPOptimizer(cfg, scorer){};

 private:
  virtual shared::optim::ILPSolver* createProblem(
      OptGraph* og, const std::set<OptNode*>& g) const;

  virtual void getConfigurationFromSolution(shared::optim::ILPSolver* lp,
                                            graph::HierarOrderCfg* c,
                                            const std::set<OptNode*>& g) const;

  void writeCrossingOracle(const std::set<OptNode*>& g,
                           shared::optim::ILPSolver* lp) const;

  void writeDiffSegConstraintsImpr(OptGraph* og, const std::set<OptNode*>& g,
                                   shared::optim::ILPSolver* lp) const;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // TRANSITMAP_OPTIM_ILPEDGEORDEROPTIMIZER_H_
