// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_ILP_ILPGRIDOPTIMIZER_H_
#define OCTI_ILP_ILPGRIDOPTIMIZER_H_

#include <vector>
#include "octi/basegraph/BaseGraph.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
#include "shared/optim/ILPSolver.h"

using octi::basegraph::BaseGraph;
using octi::basegraph::GridEdge;
using octi::basegraph::GridNode;

using octi::combgraph::CombEdge;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;

namespace octi {
namespace ilp {

struct ILPStats {
  double score;
  double time;
  size_t rows;
  size_t cols;
  bool optimal;
};

class ILPGridOptimizer {
 public:
  ILPGridOptimizer() {}

  ILPStats optimize(BaseGraph* gg, const CombGraph& cg, combgraph::Drawing* d,
                    double maxGrDist, bool noSolve,
                    const basegraph::GeoPensMap* geoPensMap, int timeLim,
                    const std::string& cacheDir, const std::string& solverStr,
                    const std::string& path) const;

 protected:
  shared::optim::ILPSolver* createProblem(
      BaseGraph* gg, const CombGraph& cg,
      const basegraph::GeoPensMap* geoPensMap, double maxGrDist,
      const std::string& solverStr) const;

  std::string getEdgUseVar(const GridEdge* e, const CombEdge* cg) const;
  std::string getStatPosVar(const GridNode* e, const CombNode* cg) const;

  void extractSolution(shared::optim::ILPSolver* lp, BaseGraph* gg,
                       const CombGraph& cg, combgraph::Drawing* d) const;

  shared::optim::StarterSol extractFeasibleSol(combgraph::Drawing* d,
                                               BaseGraph* gg,
                                               const CombGraph& cg,
                                               double maxGrDist) const;

  size_t nonInfDeg(const GridNode* g) const;
};
}  // namespace ilp
}  // namespace octi

#endif  // OCTI_ILP_ILPGRIDOPTIMIZER_H_
