// Copyright 2016
// Author: Patrick Brosi

#include <cassert>
#include <string>
#include "util/Misc.h"
#include "shared/tests/GurobiSolverTest.h"

#ifdef GUROBI_FOUND

#include "shared/optim/GurobiSolver.h"

using shared::optim::GurobiSolver;

#endif

// _____________________________________________________________________________
void GurobiSolverTest::run() {
#ifdef GUROBI_FOUND
  GurobiSolver s(shared::optim::MAX);

  size_t col1 = s.addCol("x", shared::optim::BIN, 1);
  size_t col2 = s.addCol("y", shared::optim::BIN, 1);
  size_t col3 = s.addCol("z", shared::optim::BIN, 2);

  s.update();

  TEST(s.getVarByName("x"), ==, 0);
  TEST(s.getVarByName("y"), ==, 1);
  TEST(s.getVarByName("z"), ==, 2);

  size_t row1 = s.addRow("constr1", 4, shared::optim::UP);
  s.addColToRow(col1, row1, 1);
  s.addColToRow(col2, row1, 2);
  s.addColToRow(col3, row1, 3);

  size_t row2 = s.addRow("constr2", 1, shared::optim::LO);
  s.addColToRow(col1, row2, 1);
  s.addColToRow(col2, row2, 1);

  s.update();

  TEST(s.getConstrByName("constr1"), ==, 0);
  TEST(s.getConstrByName("constr2"), ==, 1);

  s.solve();

#endif
}
