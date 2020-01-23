// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_OPTIM_GUROBISOLVER_H_
#define SHARED_OPTIM_GUROBISOLVER_H_

#ifdef GUROBI_FOUND

#include "gurobi_c.h"
#include "shared/optim/ILPSolver.h"

namespace shared {
namespace optim {

class GurobiSolver : ILPSolver {
 public:
  GurobiSolver(DirType dir);
  ~GurobiSolver();

  int addCol(const std::string& name, ColType colType, double objCoef);
  int addRow(const std::string& name, double bnd, RowType rowType);

  void addColToRow(const std::string& colName, const std::string& rowName,
                   double coef);
  void addColToRow(int colId, int rowId, double coef);

  int getVarByName(const std::string& name) const;
  int getConstrByName(const std::string& name) const;

  void solve();
  void update();

 private:
  GRBenv* _env;
  GRBmodel* _model;

  size_t _numVars, _numRows;
};

}  // namespace optim
}  // namespace shared

#endif  // SHARED_OPTIM_ILPSOLVER_H_

#endif
