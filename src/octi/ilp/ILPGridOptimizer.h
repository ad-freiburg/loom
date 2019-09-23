// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef OCTI_ILP_ILPGRIDOPTIMIZER_H_
#define OCTI_ILP_ILPGRIDOPTIMIZER_H_

#include <glpk.h>
#include <vector>
#include "octi/gridgraph/GridGraph.h"

using octi::gridgraph::GridGraph;
using octi::gridgraph::GridNode;
using octi::gridgraph::GridEdge;

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

  int optimize(GridGraph* gg) const;

 protected:
  virtual glp_prob* createProblem(const GridGraph& gg) const;

  std::string getEdgeVar(const GridEdge* e, size_t i) const;

  void solveProblem(glp_prob* lp) const;
};
}  // namespace optim
}  // namespace transitmapper

#endif  // OCTI_ILP_ILPGRIDOPTIMIZER_H_
