// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_OPTIM_ILPSOLVER_H_
#define SHARED_OPTIM_ILPSOLVER_H_

#include <string>

namespace shared {
namespace optim {

enum ColType { INT, BIN, CONT };
enum RowType { FIX, UP, LO };
enum DirType { MAX, MIN };
enum SolveType { OPTIM, INF, NON_OPTIM };

class ILPSolver {
 public:
  ILPSolver(){};

  virtual int addCol(const std::string& name, ColType colType,
                     double objCoef) = 0;
  virtual int addRow(const std::string& name, double bnd, RowType rowType) = 0;

  virtual void addColToRow(const std::string& rowName,
                           const std::string& colName, double coef) = 0;
  virtual void addColToRow(int rowId, int colId, double coef) = 0;

  virtual int getVarByName(const std::string& name) const = 0;
  virtual int getConstrByName(const std::string& name) const = 0;

  virtual void setObjCoef(const std::string& name, double coef) const = 0;
  virtual void setObjCoef(int colId, double coef) const = 0;

  virtual double getVarVal(int colId) const = 0;
  virtual double getVarVal(const std::string& name) const = 0;

  virtual SolveType solve() = 0;
  virtual void update() = 0;

  virtual double getObjVal() const = 0;

  virtual size_t getNumConstrs() const = 0;
  virtual size_t getNumVars() const = 0;
};

}  // namespace optim
}  // namespace shared

#endif  // SHARED_OPTIM_ILPSOLVER_H_
