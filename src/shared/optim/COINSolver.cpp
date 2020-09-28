// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef COIN_FOUND

#include <glpk.h>
#include <cassert>
#include <sstream>
#include <stdexcept>

// COIN includes
#include "OsiSolverInterface.hpp"
#include "OsiCbcSolverInterface.hpp"

#include "shared/optim/COINSolver.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/log/Log.h"

using shared::optim::COINSolver;
using shared::optim::SolveType;

// _____________________________________________________________________________
COINSolver::COINSolver(DirType dir)
    : _numCols(0),
      _starterArr(0),
      _status(INF),
      _timeLimit(std::numeric_limits<int>::max()) {
  LOGTO(INFO, std::cerr) << "Creating COIN solver instance...";

  assert("COIN SOLVER NOT YET IMPLEMENTED" == 0);

  // basically following
  // https://github.com/coin-or/Cbc/blob/master/examples/sample5.cpp

  _solver = &_solver1;
}

// _____________________________________________________________________________
COINSolver::~COINSolver() {
  if (_starterArr) delete[] _starterArr;
}

// _____________________________________________________________________________
int COINSolver::addCol(const std::string& name, ColType colType,
                       double objCoef) {
  return addCol(name, colType, objCoef, -COIN_DBL_MAX, COIN_DBL_MAX) + 1;
}

// _____________________________________________________________________________
int COINSolver::addCol(const std::string& name, ColType colType, double objCoef,
                       double lowBnd, double upBnd) {
  _solver->addCol(0, NULL, NULL, lowBnd, upBnd, objCoef);
  _numCols++;
  int colId = _numCols - 1;

  switch (colType) {
    case INT:
      _solver->setInteger(colId);
      break;
    case BIN:
      _solver->setInteger(colId);
      break;
    case CONT:
      _solver->setContinuous(colId);
      break;
  }

  return colId + 1;
}

// _____________________________________________________________________________
int COINSolver::addRow(const std::string& name, double bnd, RowType rowType) {
}

// _____________________________________________________________________________
void COINSolver::addColToRow(const std::string& rowName,
                             const std::string& colName, double coef) {
}

// _____________________________________________________________________________
int COINSolver::getVarByName(const std::string& name) const {
}

// _____________________________________________________________________________
int COINSolver::getConstrByName(const std::string& name) const {
}

// _____________________________________________________________________________
void COINSolver::addColToRow(int rowId, int colId, double coef) {
}

// _____________________________________________________________________________
double COINSolver::getObjVal() const {}

// _____________________________________________________________________________
SolveType COINSolver::solve() {
}

// _____________________________________________________________________________
void COINSolver::setTimeLim(int s) { _timeLimit = s * 1000; }

// _____________________________________________________________________________
int COINSolver::getTimeLim() const { return _timeLimit / 1000; }

// _____________________________________________________________________________
double* COINSolver::getStarterArr() const { return _starterArr; }

// _____________________________________________________________________________
double COINSolver::getVarVal(int colId) const {
}

// _____________________________________________________________________________
double COINSolver::getVarVal(const std::string& colName) const {
}

// _____________________________________________________________________________
void COINSolver::setObjCoef(const std::string& colName, double coef) const {
}

// _____________________________________________________________________________
void COINSolver::setObjCoef(int colId, double coef) const {
}

// _____________________________________________________________________________
void COINSolver::update() {}

// _____________________________________________________________________________
int COINSolver::getNumConstrs() const {}

// _____________________________________________________________________________
int COINSolver::getNumVars() const {}

// _____________________________________________________________________________
void COINSolver::setStarter(const StarterSol& starterSol) {
}

// _____________________________________________________________________________
void COINSolver::writeMps(const std::string& path) const {
}

#endif
