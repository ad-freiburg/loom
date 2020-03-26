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

class GurobiSolver : public ILPSolver {
 public:
  GurobiSolver(DirType dir);
  ~GurobiSolver();

  int addCol(const std::string& name, ColType colType, double objCoef);
  int addCol(const std::string& name, ColType colType, double objCoef,
             double lowBnd, double upBnd);
  int addRow(const std::string& name, double bnd, RowType rowType);

  void addColToRow(const std::string& rowName, const std::string& colName,
                   double coef);
  void addColToRow(int rowId, int colId, double coef);

  int getVarByName(const std::string& name) const;
  int getConstrByName(const std::string& name) const;

  double getVarVal(int colId) const;
  double getVarVal(const std::string& name) const;

  void setObjCoef(const std::string& name, double coef) const;
  void setObjCoef(int colId, double coef) const;

  SolveType solve();
  SolveType getStatus() { return _status; }
  void update();

  double getObjVal() const;

  int getNumConstrs() const;
  int getNumVars() const;

  void setTimeLim(int ms);
  int getTimeLim() const;

  void writeMps(const std::string& path) const;

  void setStarter(const StarterSol& starterSol);

 private:
  GRBenv* _env;
  GRBmodel* _model;

  double* _starterArr;
  StarterSol murr;

  SolveType _status;

  int _numVars, _numRows;
  std::string _logBuffer;

  static int termHook(GRBmodel* mod, void* cbdata, int where, void* solver);
};

}  // namespace optim
}  // namespace shared

#endif

#endif
