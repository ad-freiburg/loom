// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_OPTIM_COINSOLVER_H_
#define SHARED_OPTIM_COINSOLVER_H_

#ifdef COIN_FOUND

#include <vector>
#include "shared/optim/ILPSolver.h"

// COIN includes
#include "CbcSolver.hpp"
#include "CoinBuild.hpp"
#include "CoinMessageHandler.hpp"
#include "CoinModel.hpp"
#include "OsiCbcSolverInterface.hpp"
#include "OsiSolverInterface.hpp"

namespace shared {
namespace optim {

class COINSolver : public ILPSolver {
 public:
  COINSolver(DirType dir);
  ~COINSolver();

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

  void setCacheThreshold(double gb);
  double getCacheThreshold() const;

  void setNumThreads(int n);
  int getNumThreads() const;

  void setStarter(const StarterSol& starterSol);
  void writeMps(const std::string& path) const;

  double* getStarterArr() const;

 private:
  double* _starterArr;

  SolveType _status;

  int _timeLimit;

  int _numThreads;

  OsiClpSolverInterface _solver1;
  OsiSolverInterface* _solver;
  mutable CoinModel _model;
  CbcModel _cbcModel;
CoinMessageHandler _msgHandler;
};

}  // namespace optim
}  // namespace shared

#endif

#endif
