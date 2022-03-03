// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef COIN_FOUND

#include <glpk.h>
#include <cassert>
#include <sstream>
#include <stdexcept>

// COIN includes
#include "CbcSolver.hpp"
#include "CoinPragma.hpp"
#include "CoinWarmStart.hpp"
#include "OsiCbcSolverInterface.hpp"
#include "OsiSolverInterface.hpp"

#include "shared/optim/COINSolver.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/log/Log.h"

using shared::optim::COINSolver;
using shared::optim::DirType;
using shared::optim::SolveType;

// _____________________________________________________________________________
int callBack(CbcModel* model, int from) {
  int ret = 0;
  switch (from) {
    case 1:
    case 2:
      if (!model->status() && model->secondaryStatus()) ret = 1;
      break;
    case 3:
      break;
    case 4:
      break;
    case 5:
      break;
    default:
      abort();
  }
  return ret;
}

// _____________________________________________________________________________
COINSolver::COINSolver(DirType dir)
    : _starterArr(0),
      _status(INF),
      _timeLimit(std::numeric_limits<int>::max()),
      _numThreads(0),
      _msgHandler(stderr) {
  _solver = &_solver1;

  // use or own msghandler which outputs to stderr, set loglevel of CBC (=0) to
  // normal (=1)
  _msgHandler.setLogLevel(0, 1);
  _solver->passInMessageHandler(&_msgHandler);

  if (dir == MAX) {
    _model.setOptimizationDirection(-1);
    _solver->setObjSense(-1);
  } else {
    _model.setOptimizationDirection(1);
    _solver->setObjSense(1);
  }
}

// _____________________________________________________________________________
COINSolver::~COINSolver() {
  if (_starterArr) delete[] _starterArr;
}

// _____________________________________________________________________________
int COINSolver::addCol(const std::string& name, ColType colType,
                       double objCoef) {
  _model.addCol(0, NULL, NULL, -COIN_DBL_MAX, COIN_DBL_MAX, objCoef,
                name.c_str());
  int colId = _model.numberColumns() - 1;

  switch (colType) {
    case INT:
      _model.setInteger(colId);
      break;
    case BIN:
      _model.setInteger(colId);
      _model.setColLower(colId, 0.0);
      _model.setColUpper(colId, 1.0);
      break;
    case CONT:
      _model.setContinuous(colId);
      break;
  }

  return colId;
}

// _____________________________________________________________________________
int COINSolver::addCol(const std::string& name, ColType colType, double objCoef,
                       double lowBnd, double upBnd) {
  _model.addCol(0, NULL, NULL, lowBnd, upBnd, objCoef, name.c_str());
  int colId = _model.numberColumns() - 1;

  switch (colType) {
    case INT:
      _model.setInteger(colId);
      break;
    case BIN:
      _model.setInteger(colId);
      _model.setColLower(colId, 0.0);
      _model.setColUpper(colId, 1.0);
      break;
    case CONT:
      _model.setContinuous(colId);
      break;
  }

  return colId;
}

// _____________________________________________________________________________
int COINSolver::addRow(const std::string& name, double bnd, RowType rowType) {
  switch (rowType) {
    case FIX:
      _model.addRow(0, 0, 0, bnd, bnd, name.c_str());
      break;
    case UP:
      _model.addRow(0, 0, 0, -COIN_DBL_MAX, bnd, name.c_str());
      break;
    case LO:
      _model.addRow(0, 0, 0, bnd, COIN_DBL_MAX, name.c_str());
      break;
  }

  int rowId = _model.numberRows() - 1;

  return rowId;
}

// _____________________________________________________________________________
void COINSolver::addColToRow(const std::string& rowName,
                             const std::string& colName, double coef) {
  addColToRow(getConstrByName(rowName), getVarByName(colName), coef);
}

// _____________________________________________________________________________
int COINSolver::getVarByName(const std::string& name) const {
  return _model.column(name.c_str());
}

// _____________________________________________________________________________
int COINSolver::getConstrByName(const std::string& name) const {
  return _model.row(name.c_str());
}

// _____________________________________________________________________________
void COINSolver::addColToRow(int rowId, int colId, double coef) {
  _model.setElement(rowId, colId, coef);
}

// _____________________________________________________________________________
double COINSolver::getObjVal() const { return _solver->getObjValue(); }

// _____________________________________________________________________________
SolveType COINSolver::solve() {
  _solver->loadFromCoinModel(_model);

  _solver1.getModelPtr()->setMoreSpecialOptions(3);
  _cbcModel = CbcModel(_solver1);

  _cbcModel.setMaximumSeconds(_timeLimit);
  _cbcModel.setUseElapsedTime(true);

  // this basically follows the examle given at
  // https://github.com/coin-or/Cbc/blob/879602724a65987c5cfb0b9fbacfa192c93df42c/examples/driver6.cpp
  //
  // It's a bit finicky to use the CBC solver. If the library functions are
  // called directly, the user has to set up the solver with sensible default
  // values, add all heuristics etc. with sensible default values, and start the
  // solving process. This is different than in other libraries (gurobi, glpk),
  // where the solver offers methods which uses carefully tuned default settings
  // (the settings which are used when the command line interface is called
  // directly). To get this behavior with CBC, we simulate the process used by
  // the cbc command line interface and call CbcMain0 and CbcMain1. CbcMain1
  // expects exactly the same arguments (in array argv) as the command line
  // interface. The code below therefor acts like if we call "cbc -solve
  // -threads <N> ilp.mps" from the command line (verified by several tests on
  // large ILPs which yield basically the same solution time and command line
  // output).

  CbcSolverUsefulData solverData;
  CbcMain0(_cbcModel, solverData);
  std::string numThreads = "4";

  if (_numThreads > 0) numThreads = std::to_string(_numThreads);

  const char* argv2[] = {"-solve", "-threads", numThreads.c_str()};
  CbcMain1(3, argv2, _cbcModel, callBack, solverData);
  _solver = _cbcModel.solver();

  if (_cbcModel.isProvenOptimal())
    _status = OPTIM;
  else if (_cbcModel.isProvenInfeasible())
    _status = INF;
  else if (_cbcModel.bestSolution())
    _status = NON_OPTIM;
  else
    _status = INF;

  return getStatus();
}

// _____________________________________________________________________________
void COINSolver::setTimeLim(int s) { _timeLimit = s; }

// _____________________________________________________________________________
int COINSolver::getTimeLim() const { return _timeLimit; }

// _____________________________________________________________________________
double* COINSolver::getStarterArr() const { return _starterArr; }

// _____________________________________________________________________________
double COINSolver::getVarVal(int colId) const {
  return _solver->getColSolution()[colId];
}

// _____________________________________________________________________________
double COINSolver::getVarVal(const std::string& colName) const {
  return getVarVal(getVarByName(colName));
}

// _____________________________________________________________________________
void COINSolver::setObjCoef(const std::string& colName, double coef) const {
  setObjCoef(getVarByName(colName), coef);
}

// _____________________________________________________________________________
void COINSolver::setObjCoef(int colId, double coef) const {
  _model.setObjective(colId, coef);
}

// _____________________________________________________________________________
void COINSolver::update() {}

// _____________________________________________________________________________
int COINSolver::getNumConstrs() const { return _model.numberRows(); }

// _____________________________________________________________________________
int COINSolver::getNumVars() const { return _model.numberColumns(); }

// _____________________________________________________________________________
void COINSolver::setStarter(const StarterSol& starterSol) {
  // TODO: not yet implemented
  UNUSED(starterSol);
}

// _____________________________________________________________________________
void COINSolver::setNumThreads(int n) {
  LOGTO(INFO, std::cerr) << "Setting number of threads to " << n;
  _numThreads = n;
}

// _____________________________________________________________________________
int COINSolver::getNumThreads() const { return _numThreads; }

// _____________________________________________________________________________
void COINSolver::writeMps(const std::string& path) const {
  _model.writeMps(path.c_str());
}

// _____________________________________________________________________________
void COINSolver::setCacheDir(const std::string& dir) {
  LOGTO(INFO, std::cerr) << "Setting cache dir to " << dir
                         << " (TODO: not implemented for COIN)";
}

// _____________________________________________________________________________
void COINSolver::setCacheThreshold(double gb) {
  LOGTO(INFO, std::cerr) << "Setting cache threshhold to " << gb
                         << " GB (TODO: not implemented for COIN)";
}

// _____________________________________________________________________________
double COINSolver::getCacheThreshold() const { return 0; }

// _____________________________________________________________________________
std::string COINSolver::getCacheDir() const { return ""; }

#endif
