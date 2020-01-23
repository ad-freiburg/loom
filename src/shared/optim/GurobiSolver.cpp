// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef GUROBI_FOUND

#include "gurobi_c.h"
#include "shared/optim/GurobiSolver.h"
#include "util/log/Log.h"

using shared::optim::GurobiSolver;

// _____________________________________________________________________________
GurobiSolver::GurobiSolver(DirType dir) : _numVars(0), _numRows(0) {
  LOG(DEBUG) << "Creating gurobi solver instance...";

  int error = GRBloadenv(&_env, "gurobi.log");
  if (error || _env == 0) {
    LOG(ERROR) << "Could not create gurobi environment";
    exit(1);
  }

  // create emtpy model
  error = GRBnewmodel(_env, &_model, "mip1", 0, 0, 0, 0, 0, 0);
  if (error) {
    LOG(ERROR) << "Could not create gurobi model";
    exit(1);
  }

  if (dir == MAX)
    GRBsetintattr(_model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
  else
    GRBsetintattr(_model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
}

// _____________________________________________________________________________
GurobiSolver::~GurobiSolver() {
  LOG(DEBUG) << "Destroying gurobi solver...";
  GRBfreeenv(_env);
}

// _____________________________________________________________________________
int GurobiSolver::addCol(const std::string& name, ColType colType,
                            double objCoef) {
  char vtype;
  switch (colType) {
    case INT:
      vtype = GRB_INTEGER;
      break;
    case BIN:
      vtype = GRB_BINARY;
      break;
    case CONT:
      vtype = GRB_CONTINUOUS;
      break;
  }
  int error = GRBaddvar(_model, 0, 0, 0, objCoef, -GRB_INFINITY, GRB_INFINITY,
                        vtype, name.c_str());
  if (error) {
    LOG(ERROR) << "Could not add variable " << name;
    exit(1);
  }
  _numVars++;
  return _numVars - 1;
}

// _____________________________________________________________________________
int GurobiSolver::addRow(const std::string& name, double bnd,
                            RowType rowType) {
  char rtype;
  switch (rowType) {
    case FIX:
      rtype = GRB_EQUAL;
      break;
    case UP:
      rtype = GRB_LESS_EQUAL;
      break;
    case LO:
      rtype = GRB_GREATER_EQUAL;
      break;
  }
  int error = GRBaddconstr(_model, 0, 0, 0, rtype, bnd, name.c_str());
  if (error) {
    LOG(ERROR) << "Could not add row " << name;
    exit(1);
  }

  _numRows++;
  return _numRows - 1;
}

// _____________________________________________________________________________
void GurobiSolver::addColToRow(const std::string& colName,
                               const std::string& rowName, double coef) {}

// _____________________________________________________________________________
int GurobiSolver::getVarByName(const std::string& name) const {
  int ret;
  int error = GRBgetvarbyname(_model, name.c_str(), &ret);
  if (error) return -1;
  return ret;
}

// _____________________________________________________________________________
int GurobiSolver::getConstrByName(const std::string& name) const {
  int ret;
  int error = GRBgetconstrbyname(_model, name.c_str(), &ret);
  if (error) return -1;
  return ret;
}

// _____________________________________________________________________________
void GurobiSolver::addColToRow(int colId, int rowId, double coef) {
  int col = colId;
  int row = rowId;

  int error = GRBchgcoeffs(_model, 1, &row, &col, &coef);
  if (error) {
    LOG(ERROR) << "Could not add col " << colId << " to row " << rowId;
    exit(1);
  }
}

// _____________________________________________________________________________
void GurobiSolver::solve() {
  update();
  GRBwrite(_model, "mip1.lp");
  int error = GRBoptimize(_model);
  if (error) {
    LOG(ERROR) << "Could not optimize model";
    exit(1);
  }

  int status;
  error = GRBgetintattr(_model, GRB_INT_ATTR_STATUS, &status);
  if (error) {
    LOG(ERROR) << "Could not retrieve optim status";
    exit(1);
  }
}

// _____________________________________________________________________________
void GurobiSolver::update() { GRBupdatemodel(_model); }

#endif
