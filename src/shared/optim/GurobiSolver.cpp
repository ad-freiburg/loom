// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef GUROBI_FOUND

#include <sstream>
#include <stdexcept>
#include "gurobi_c.h"
#include "shared/optim/GurobiSolver.h"
#include "util/log/Log.h"

using shared::optim::GurobiSolver;
using shared::optim::SolveType;

// _____________________________________________________________________________
GurobiSolver::GurobiSolver(DirType dir)
    : _status(INF), _numVars(0), _numRows(0) {
  LOG(DEBUG) << "Creating gurobi solver instance...";

  int error = GRBloadenv(&_env, "gurobi.log");
  if (error || _env == 0) {
    throw std::runtime_error("Could not create gurobi environment");
  }


  // create emtpy model
  error = GRBnewmodel(_env, &_model, "loom_mip", 0, 0, 0, 0, 0, 0);
  if (error) {
    throw std::runtime_error("Could not create gurobi model");
  }

  if (dir == MAX)
    GRBsetintattr(_model, GRB_INT_ATTR_MODELSENSE, GRB_MAXIMIZE);
  else
    GRBsetintattr(_model, GRB_INT_ATTR_MODELSENSE, GRB_MINIMIZE);
}

// _____________________________________________________________________________
GurobiSolver::~GurobiSolver() {
  GRBfreemodel(_model);
  GRBfreeenv(_env);
  if (_starterArr) delete[] _starterArr;
}

// _____________________________________________________________________________
int GurobiSolver::addCol(const std::string& name, ColType colType,
                         double objCoef) {
  return addCol(name, colType, objCoef, -GRB_INFINITY, GRB_INFINITY);
}

// _____________________________________________________________________________
int GurobiSolver::addCol(const std::string& name, ColType colType,
                         double objCoef, double lowBnd, double upBnd) {
  char vtype = 0;
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
  int error =
      GRBaddvar(_model, 0, 0, 0, objCoef, lowBnd, upBnd, vtype, name.c_str());
  if (error) {
    throw std::runtime_error("Could not add variable " + name);
  }
  _numVars++;
  return _numVars - 1;
}

// _____________________________________________________________________________
int GurobiSolver::addRow(const std::string& name, double bnd, RowType rowType) {
  char rtype = 0;
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
    throw std::runtime_error("Could not add row " + name);
  }

  _numRows++;
  return _numRows - 1;
}

// _____________________________________________________________________________
void GurobiSolver::addColToRow(const std::string& rowName,
                               const std::string& colName, double coef) {
  int col = getVarByName(colName);
  if (col < 0) {
    LOG(ERROR) << "Could not find variable " << colName;
  }
  int row = getConstrByName(rowName);
  if (row < 0) {
    LOG(ERROR) << "Could not find constraint " << rowName;
  }

  addColToRow(col, row, coef);
}

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
void GurobiSolver::addColToRow(int rowId, int colId, double coef) {
  int col = colId;
  int row = rowId;

  int error = GRBchgcoeffs(_model, 1, &row, &col, &coef);
  if (error) {
    std::stringstream ss;
    ss << "Could not add col " << colId << " to row " << rowId << " (" << error
       << ")";
    throw std::runtime_error(ss.str());
  }
}

// _____________________________________________________________________________
double GurobiSolver::getObjVal() const {
  double objVal;

  int error = GRBgetdblattr(_model, GRB_DBL_ATTR_OBJVAL, &objVal);
  if (error) {
    throw std::runtime_error(
        "Could not retrieve optimal target function value.");
  }

  return objVal;
}

// _____________________________________________________________________________
void GurobiSolver::setStarter(const StarterSol& starterSol) {
  _starterArr = new double[getNumVars()];
  std::fill_n(_starterArr, getNumVars(), GRB_UNDEFINED);

  for (const auto& varVal : starterSol) {
    int colId = getVarByName(varVal.first);
    if (colId < 0) continue;
    _starterArr[colId] = varVal.second;
  }
}

// _____________________________________________________________________________
SolveType GurobiSolver::solve() {
  update();

  int error;

  // error = GRBtunemodel(_model);
  // int nresults;

  // GRBsetdblparam(_env, "TuneTimeLimit", 60.0 * 1);

  // error = GRBgetintattr(_model, "TuneResultCount", &nresults);
  // if (nresults > 0) {
  // error = GRBgettuneresult(_model, 0);
  // }

  if (_starterArr) {
    error = GRBsetdblattrarray(_model, "Start", 0, getNumVars(), _starterArr);
    if (error) {
      throw std::runtime_error("Could not set start solution");
    }
  }

  error = GRBoptimize(_model);
  if (error) {
    throw std::runtime_error("Could not optimize model");
  }

  int status;
  error = GRBgetintattr(_model, GRB_INT_ATTR_STATUS, &status);
  if (error) {
    throw std::runtime_error("Could not retrieve optim status");
  }

  int optimStat;
  double objVal;

  error = GRBgetintattr(_model, GRB_INT_ATTR_STATUS, &optimStat);
  if (error) {
    throw std::runtime_error("Could not read optimization status.");
  }

  error = GRBgetdblattr(_model, GRB_DBL_ATTR_OBJVAL, &objVal);
  if (error) {
    throw std::runtime_error(
        "Could not retrieve optimal target function value.");
  }

  if (optimStat == GRB_OPTIMAL) _status = OPTIM;
  else if (optimStat == GRB_INF_OR_UNBD) _status = INF;
  else _status = NON_OPTIM;

  return getStatus();
}

// _____________________________________________________________________________
void GurobiSolver::setTimeLim(int s) {
  // set time limit
  GRBsetdblparam(_env, "TimeLimit", s);
}

// _____________________________________________________________________________
int GurobiSolver::getTimeLim() const {
  double ret;
  int error = GRBgetdblparam(_env, "TimeLimit", &ret);
  if (error) {
    std::stringstream ss;
    ss << "Could not retrieve time limit value";
    throw std::runtime_error(ss.str());
  }
  return ret;
}

// _____________________________________________________________________________
double GurobiSolver::getVarVal(int colId) const {
  double val;
  int error = GRBgetdblattrelement(_model, GRB_DBL_ATTR_X, colId, &val);
  if (error) {
    std::stringstream ss;
    ss << "Could not retrieve value for field " << colId;
    throw std::runtime_error(ss.str());
  }
  return val;
}

// _____________________________________________________________________________
double GurobiSolver::getVarVal(const std::string& colName) const {
  int col = getVarByName(colName);
  if (col < 0) {
    LOG(ERROR) << "Could not find variable " << colName;
  }

  return getVarVal(col);
}

// _____________________________________________________________________________
void GurobiSolver::setObjCoef(const std::string& colName, double coef) const {
  int col = getVarByName(colName);
  if (col < 0) {
    LOG(ERROR) << "Could not find variable " << colName;
  }

  setObjCoef(col, coef);
}

// _____________________________________________________________________________
void GurobiSolver::setObjCoef(int colId, double coef) const {
  int error = GRBsetdblattrelement(_model, GRB_DBL_ATTR_OBJ, colId, coef);
  if (error) {
    std::stringstream ss;
    ss << "Could not change objective value for col " << colId;
    throw std::runtime_error(ss.str());
  }
}

// _____________________________________________________________________________
void GurobiSolver::update() { GRBupdatemodel(_model); }

// _____________________________________________________________________________
int GurobiSolver::getNumConstrs() const { return _numRows; }

// _____________________________________________________________________________
int GurobiSolver::getNumVars() const { return _numVars; }

// _____________________________________________________________________________
void GurobiSolver::writeMps(const std::string& path) const {
  // TODO: exception if could not be written
  GRBwrite(_model, path.c_str());
}

#endif
