// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef GUROBI_FOUND

#include "gurobi_c.h"
#include "shared/optim/GurobiSolver.h"
#include "util/log/Log.h"

using shared::optim::GurobiSolver;

// _____________________________________________________________________________
GurobiSolver::GurobiSolver() {
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
}

// _____________________________________________________________________________
GurobiSolver::~GurobiSolver() {
  LOG(DEBUG) << "Destroying gurobi solver...";
  GRBfreeenv(_env);
}

// _____________________________________________________________________________
size_t GurobiSolver::addCol(const std::string& name, ColType colType,
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
  GRBaddvar(_model, 0, 0, 0, objCoef, 0, 0, vtype, name.c_str());
  _numVars++;
  return _numVars - 1;
}

// _____________________________________________________________________________
size_t GurobiSolver::addRow(const std::string& name, RowType rowType) {

}

// _____________________________________________________________________________
void GurobiSolver::addColToRow(const std::string& colName,
                               const std::string& rowName, double coef) {}

// _____________________________________________________________________________
void GurobiSolver::addColToRow(size_t colId, size_t rowId, double coef) {
  // GRBchgcoef does the right thing here
}

// _____________________________________________________________________________
void GurobiSolver::solve() {}

// _____________________________________________________________________________
void GurobiSolver::update() {}

#endif
