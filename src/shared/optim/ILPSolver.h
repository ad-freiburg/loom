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

class ILPSolver {
 public:
  ILPSolver(){};

  virtual size_t addCol(const std::string& name, ColType colType,
                        double objCoef) = 0;
  virtual size_t addRow(const std::string& name, RowType rowType) = 0;

  virtual void addColToRow(const std::string& colName,
                           const std::string& rowName, double coef) = 0;
  virtual void addColToRow(size_t colId, size_t rowId, double coef) = 0;

  virtual void solve() = 0;
  virtual void update() = 0;
};

}  // namespace optim
}  // namespace shared

#endif  // SHARED_OPTIM_ILPSOLVER_H_
