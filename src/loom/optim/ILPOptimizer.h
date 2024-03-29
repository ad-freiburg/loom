// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef LOOM_OPTIM_ILPOPTIMIZER_H_
#define LOOM_OPTIM_ILPOPTIMIZER_H_

#include "loom/optim/ExhaustiveOptimizer.h"
#include "loom/config/LoomConfig.h"
#include "loom/optim/OptGraph.h"
#include "loom/optim/Optimizer.h"
#include "shared/linegraph/Line.h"
#include "shared/optim/ILPSolver.h"
#include "shared/rendergraph/OrderCfg.h"

namespace loom {
namespace optim {

class ILPOptimizer : public Optimizer {
 public:
  ILPOptimizer(const config::Config* cfg,
               const shared::rendergraph::Penalties& pens)
      : Optimizer(cfg, pens), _exhausOpt(cfg, pens) {};

  virtual double optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                              shared::rendergraph::HierarOrderCfg* c,
                              size_t depth, OptResStats& stats) const;

  virtual std::string getName() const { return "ilp";}

 protected:
  const loom::optim::ExhaustiveOptimizer _exhausOpt;
  virtual shared::optim::ILPSolver* createProblem(
      OptGraph* og, const std::set<OptNode*>& g) const;

  virtual void getConfigurationFromSolution(
      shared::optim::ILPSolver* lp, shared::rendergraph::HierarOrderCfg* c,
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
  int getSeparationPenalty(const OptNode* n) const;

  bool separationOpt() const;
};
}  // namespace optim
}  // namespace loom

#endif  // LOOM_OPTIM_ILPOPTIMIZER_H_
