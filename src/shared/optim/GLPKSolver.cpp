// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef GLPK_FOUND

#include <glpk.h>
#include <cassert>
#include <sstream>
#include <stdexcept>
#include "shared/optim/GLPKSolver.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/log/Log.h"

using shared::optim::GLPKSolver;
using shared::optim::SolveType;
using shared::optim::VariableMatrix;

// _____________________________________________________________________________
GLPKSolver::GLPKSolver(DirType dir)
    : _starterArr(0),
      _status(INF),
      _timeLimit(std::numeric_limits<int>::max()) {
  LOG(INFO, std::cerr) << "Creating GLPK solver instance...";

  _prob = glp_create_prob();

  glp_set_prob_name(_prob, "loom_mip");

  if (dir == MAX)
    glp_set_obj_dir(_prob, GLP_MAX);
  else
    glp_set_obj_dir(_prob, GLP_MIN);

  glp_create_index(_prob);
}

// _____________________________________________________________________________
GLPKSolver::~GLPKSolver() {
  glp_delete_prob(_prob);
  glp_free_env();
  if (_starterArr) delete[] _starterArr;
}

// _____________________________________________________________________________
int GLPKSolver::addCol(const std::string& name, ColType colType,
                       double objCoef) {
  int vtype = 0;
  switch (colType) {
    case INT:
      vtype = GLP_IV;
      break;
    case BIN:
      vtype = GLP_BV;
      break;
    case CONT:
      vtype = GLP_CV;
      break;
  }

  int col = glp_add_cols(_prob, 1);
  glp_set_col_name(_prob, col, name.c_str());
  glp_set_col_kind(_prob, col, vtype);
  glp_set_obj_coef(_prob, col, objCoef);

  return col - 1;
}

// _____________________________________________________________________________
int GLPKSolver::addCol(const std::string& name, ColType colType, double objCoef,
                       double lowBnd, double upBnd) {
  int rtype = 0;
  if (lowBnd <= -std::numeric_limits<double>::max() &&
      upBnd >= std::numeric_limits<double>::max()) {
    rtype = GLP_FR;
  } else if (lowBnd <= -std::numeric_limits<double>::max()) {
    rtype = GLP_UP;
  } else if (upBnd >= std::numeric_limits<double>::max()) {
    rtype = GLP_LO;
  } else if (lowBnd == upBnd) {
    rtype = GLP_FX;
  } else {
    rtype = GLP_DB;
  }

  int col = addCol(name, colType, objCoef);
  glp_set_col_bnds(_prob, col + 1, rtype, lowBnd, upBnd);

  return col;
}

// _____________________________________________________________________________
int GLPKSolver::addRow(const std::string& name, double bnd, RowType rowType) {
  int rtype = 0;
  switch (rowType) {
    case FIX:
      rtype = GLP_FX;
      break;
    case UP:
      rtype = GLP_UP;
      break;
    case LO:
      rtype = GLP_LO;
      break;
  }

  int row = glp_add_rows(_prob, 1);
  assert(row);
  glp_set_row_name(_prob, row, name.c_str());
  glp_set_row_bnds(_prob, row, rtype, bnd, bnd);

  return row - 1;
}

// _____________________________________________________________________________
void GLPKSolver::addColToRow(const std::string& rowName,
                             const std::string& colName, double coef) {
  int col = getVarByName(colName);
  if (col < 0) {
    LOG(ERROR, std::cerr) << "Could not find variable " << colName;
  }
  int row = getConstrByName(rowName);
  if (row < 0) {
    LOG(ERROR, std::cerr) << "Could not find constraint " << rowName;
  }

  addColToRow(col, row, coef);
}

// _____________________________________________________________________________
int GLPKSolver::getVarByName(const std::string& name) const {
  int ret = glp_find_col(_prob, name.c_str());
  if (!ret) return -1;
  return ret - 1;
}

// _____________________________________________________________________________
int GLPKSolver::getConstrByName(const std::string& name) const {
  int ret = glp_find_row(_prob, name.c_str());
  if (!ret) return -1;
  return ret - 1;
}

// _____________________________________________________________________________
void GLPKSolver::addColToRow(int rowId, int colId, double coef) {
  _vm.addVar(rowId + 1, colId + 1, coef);
}

// _____________________________________________________________________________
double GLPKSolver::getObjVal() const { return glp_mip_obj_val(_prob); }

// _____________________________________________________________________________
SolveType GLPKSolver::solve() {
  update();

  int* ia = 0;
  int* ja = 0;
  double* res = 0;
  _vm.getGLPKArrs(&ia, &ja, &res);

  glp_load_matrix(_prob, _vm.getNumVars(), ia, ja, res);

  delete[](ia);
  delete[](ja);
  delete[](res);

  std::string buf;
  glp_term_hook(termHook, &buf);

  glp_iocp params;
  glp_smcp sparams;

  // default initialization
  glp_init_iocp(&params);
  glp_init_smcp(&sparams);

  params.cb_func = optCb;
  params.cb_info = this;

  // params.presolve = GLP_ON;
  // params.binarize = GLP_OFF;
  // params.ps_tm_lim = 10000;
  params.tm_lim = _timeLimit;
  // params.fp_heur = GLP_ON;
  // params.ps_heur = GLP_ON;

  glp_simplex(_prob, &sparams);
  glp_intopt(_prob, &params);

  int optimStat = glp_get_status(_prob);

  if (optimStat == GLP_OPT) _status = OPTIM;
  else if (optimStat == GLP_NOFEAS || optimStat == GLP_INFEAS ||
      optimStat == GLP_UNBND || optimStat == GLP_UNDEF)
    _status = INF;
  else _status = NON_OPTIM;

  return getStatus();
}

// _____________________________________________________________________________
void GLPKSolver::setTimeLim(int s) { _timeLimit = s * 1000; }

// _____________________________________________________________________________
int GLPKSolver::getTimeLim() const { return _timeLimit / 1000; }

// _____________________________________________________________________________
double* GLPKSolver::getStarterArr() const { return _starterArr; }

// _____________________________________________________________________________
void GLPKSolver::optCb(glp_tree* tree, void* solver) {
  auto _this = reinterpret_cast<GLPKSolver*>(solver);
  switch (glp_ios_reason(tree)) {
    case GLP_IHEUR:
      if (_this->getStarterArr()) {
        glp_ios_heur_sol(tree, _this->getStarterArr());
      }
      break;
    default:
      break;
  }
}

// _____________________________________________________________________________
int GLPKSolver::termHook(void* info, const char* str) {
  std::string* buff = reinterpret_cast<std::string*>(info);
  std::string s = str;
  for (auto ch : s) {
    if (ch == '\n') {
      std::string trimmed = util::ltrim(*buff);
      if (trimmed.size()) LOG(INFO, std::cerr) << trimmed;
      buff->clear();
    } else {
      buff->push_back(ch);
    }
  }
  return 1;
}

// _____________________________________________________________________________
double GLPKSolver::getVarVal(int colId) const {
  double val = glp_mip_col_val(_prob, colId + 1);
  return val;
}

// _____________________________________________________________________________
double GLPKSolver::getVarVal(const std::string& colName) const {
  int col = getVarByName(colName);
  if (col < 0) {
    LOG(ERROR, std::cerr) << "Could not find variable " << colName;
  }

  return getVarVal(col);
}

// _____________________________________________________________________________
void GLPKSolver::setObjCoef(const std::string& colName, double coef) const {
  int col = getVarByName(colName);
  if (col < 0) {
    LOG(ERROR, std::cerr) << "Could not find variable " << colName;
  }

  setObjCoef(col, coef);
}

// _____________________________________________________________________________
void GLPKSolver::setObjCoef(int colId, double coef) const {
  glp_set_obj_coef(_prob, colId + 1, coef);
}

// _____________________________________________________________________________
void GLPKSolver::update() { glp_create_index(_prob); }

// _____________________________________________________________________________
int GLPKSolver::getNumConstrs() const { return glp_get_num_rows(_prob); }

// _____________________________________________________________________________
int GLPKSolver::getNumVars() const { return glp_get_num_cols(_prob); }

// _____________________________________________________________________________
void GLPKSolver::setStarter(const StarterSol& starterSol) {
  _starterArr = new double[getNumVars() + 1];

  size_t a = 0;

  for (const auto& varVal : starterSol) {
    int colId = getVarByName(varVal.first);
    if (colId < 0) continue;
    _starterArr[colId + 1] = varVal.second;
    a++;
  }
}

// _____________________________________________________________________________
void VariableMatrix::addVar(int row, int col, double val) {
  rowNum.push_back(row);
  colNum.push_back(col);
  vals.push_back(val);
}

// _____________________________________________________________________________
void VariableMatrix::getGLPKArrs(int** ia, int** ja, double** r) const {
  assert(rowNum.size() == colNum.size());
  assert(colNum.size() == vals.size());

  *ia = new int[rowNum.size() + 1];
  *ja = new int[rowNum.size() + 1];
  *r = new double[rowNum.size() + 1];

  // glpk arrays always start at 1 for some reason
  for (size_t i = 1; i <= rowNum.size(); ++i) {
    (*ia)[i] = rowNum[i - 1];
    (*ja)[i] = colNum[i - 1];
    (*r)[i] = vals[i - 1];
  }
}

// _____________________________________________________________________________
void GLPKSolver::writeMps(const std::string& path) const {
  // TODO: exception if could not be written
  std::string buf;
  glp_term_hook(termHook, &buf);
  glp_write_mps(_prob, GLP_MPS_FILE, 0, path.c_str());
}

#endif
