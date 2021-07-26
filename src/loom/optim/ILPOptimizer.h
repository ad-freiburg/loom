// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_ILPOPTIMIZER_H_
#define LOOM_OPTIM_ILPOPTIMIZER_H_

#include "shared/linegraph/Line.h"
#include "shared/optim/ILPSolver.h"
#include "loom/config/TransitMapConfig.h"
#include "loom/graph/OrderCfg.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/Optimizer.h"
#include "loom/optim/Scorer.h"

namespace loom {
namespace optim {

using namespace graph;

class ILPOptimizer : public Optimizer {
 public:
  ILPOptimizer(const config::Config* cfg, const Scorer* scorer)
      : Optimizer(cfg, scorer){};

  int optimizeComp(OptGraph* og, const std::set<OptNode*>& g, HierarOrderCfg* c,
                   size_t depth) const;

 protected:
  virtual shared::optim::ILPSolver* createProblem(
      OptGraph* og, const std::set<OptNode*>& g) const;

  virtual void getConfigurationFromSolution(shared::optim::ILPSolver* lp,
                                            HierarOrderCfg* c,
                                            const std::set<OptNode*>& g) const;

  std::string getILPVarName(OptEdge* e, const shared::linegraph::Line* r,
                            size_t p) const;

  void writeSameSegConstraints(OptGraph* og, const std::set<OptNode*>& g,
                               shared::optim::ILPSolver* lp) const;

  void writeDiffSegConstraints(OptGraph* og, const std::set<OptNode*>& g,
                               shared::optim::ILPSolver* lp) const;

  std::vector<PosComPair> getPositionCombinations(OptEdge* a, OptEdge* b) const;
  std::vector<PosCom> getPositionCombinations(OptEdge* a) const;

  int getCrossingPenaltySameSeg(const OptNode* n) const;
  int getCrossingPenaltyDiffSeg(const OptNode* n) const;
  int getSplittingPenalty(const OptNode* n) const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_ILPOPTIMIZER_H_
