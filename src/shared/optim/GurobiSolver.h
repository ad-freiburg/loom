// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_OPTIM_GUROBISOLVER_H_
#define SHARED_OPTIM_GUROBISOLVER_H_

#include "gurobi_c.h"
#include "shared/optim/ILPSolver.h"

namespace shared {
namespace optim {

class GurobiSolver : ILPSolver {
 public:
  GurobiSolver();
  ~GurobiSolver();

  size_t addCol(const std::string& name, ColType colType,
                        double objCoef);
  size_t addRow(const std::string& name, RowType rowType);

  void addColToRow(const std::string& colName,
                           const std::string& rowName, double coef);
  void addColToRow(size_t colId, size_t rowId, double coef);

  void solve();
  void update();

 private:
  GRBenv* _env;
  GRBmodel* _model;

  size_t _numVars;
};

}  // namespace optim
}  // namespace shared

#endif  // SHARED_OPTIM_ILPSOLVER_H_
