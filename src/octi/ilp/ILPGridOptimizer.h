// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_ILP_ILPGRIDOPTIMIZER_H_
#define OCTI_ILP_ILPGRIDOPTIMIZER_H_

#include <vector>
#include "octi/combgraph/CombGraph.h"
#include "octi/combgraph/Drawing.h"
#include "octi/gridgraph/GridGraph.h"
#include "shared/optim/ILPSolver.h"

using octi::gridgraph::GridEdge;
using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;

using octi::combgraph::CombEdge;
using octi::combgraph::CombGraph;
using octi::combgraph::CombNode;

namespace octi {
namespace ilp {

struct VariableMatrix {
  std::vector<int> rowNum;
  std::vector<int> colNum;
  std::vector<double> vals;

  void addVar(int row, int col, double val);
  void getGLPKArrs(int** ia, int** ja, double** r) const;
  size_t getNumVars() const { return vals.size(); }
};

class ILPGridOptimizer {
 public:
  ILPGridOptimizer() {}

  double optimize(GridGraph* gg, const CombGraph& cg, combgraph::Drawing* d,
                  double maxGrDist, bool noSolve,
                  const gridgraph::GeoPensMap* geoPensMap,
                  const std::string& path) const;

 protected:
  shared::optim::ILPSolver* createProblem(
      GridGraph* gg, const CombGraph& cg,
      const gridgraph::GeoPensMap* geoPensMap, double maxGrDist) const;

  std::string getEdgeUseVar(const GridEdge* e, const CombEdge* cg) const;
  std::string getStatPosVar(const GridNode* e, const CombNode* cg) const;

  void extractSolution(shared::optim::ILPSolver* lp, GridGraph* gg,
                       const CombGraph& cg, combgraph::Drawing* d) const;

  shared::optim::StarterSol extractFeasibleSol(GridGraph* gg,
                                               const CombGraph& cg,
                                               double maxGrDist) const;

  size_t nonInfDeg(const GridNode* g) const;
};
}  // namespace ilp
}  // namespace octi

#endif  // OCTI_ILP_ILPGRIDOPTIMIZER_H_
