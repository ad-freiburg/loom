// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_OPTIM_GLPKSOLVER_H_
#define SHARED_OPTIM_GLPKSOLVER_H_

#ifdef GLPK_FOUND

#include <glpk.h>
#include <vector>
#include "shared/optim/ILPSolver.h"

namespace shared {
namespace optim {

struct VariableMatrix {
  std::vector<int> rowNum;
  std::vector<int> colNum;
  std::vector<double> vals;

  void addVar(int row, int col, double val);
  void getGLPKArrs(int** ia, int** ja, double** r) const;
  size_t getNumVars() const { return vals.size(); }
};

class GLPKSolver : public ILPSolver {
 public:
  GLPKSolver(DirType dir);
  ~GLPKSolver();

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

  void setTimeLim(int s);
  int getTimeLim() const;

  void setCacheDir(const std::string& dir);
  std::string getCacheDir() const;

  void setStarter(const StarterSol& starterSol);
  void writeMps(const std::string& path) const;

  double* getStarterArr() const;

 private:
  glp_prob* _prob;
  VariableMatrix _vm;

  double* _starterArr;

  SolveType _status;

  int _timeLimit;

  std::string _termBuf;

  static void optCb(glp_tree* tree, void* solver);
  static int termHook(void* info, const char* str);
  static void errorHook(void* info);
};

}  // namespace optim
}  // namespace shared

#endif

#endif
