// Copyright 2016
// Author: Patrick Brosi

#include <cassert>
#include <string>
#include "shared/tests/GurobiSolverTest.h"
#include "util/Misc.h"

#ifdef GUROBI_FOUND

#include "shared/optim/GurobiSolver.h"

using shared::optim::GurobiSolver;
using util::approx;

#endif

// _____________________________________________________________________________
void GurobiSolverTest::run() {
#ifdef GUROBI_FOUND

  {
    GurobiSolver s(shared::optim::MAX);

    size_t col1 = s.addCol("x", shared::optim::BIN, 1);
    size_t col2 = s.addCol("y", shared::optim::BIN, 1);
    size_t col3 = s.addCol("z", shared::optim::BIN, 2);

    s.update();

    TEST(s.getVarByName("x"), ==, 0);
    TEST(s.getVarByName("y"), ==, 1);
    TEST(s.getVarByName("z"), ==, 2);

    size_t row1 = s.addRow("constr1", 4, shared::optim::UP);
    s.addColToRow(row1, col1, 1);
    s.addColToRow(row1, col2, 2);
    s.addColToRow(row1, col3, 3);

    size_t row2 = s.addRow("constr2", 1, shared::optim::LO);
    s.addColToRow(row2, col1, 1);
    s.addColToRow(row2, col2, 1);

    s.update();

    TEST(s.getConstrByName("constr1"), ==, 0);
    TEST(s.getConstrByName("constr2"), ==, 1);

    auto ret = s.solve();

    TEST(ret, ==, shared::optim::OPTIM);

    TEST(s.getVarVal("x"), ==, approx(1));
    TEST(s.getVarVal("y"), ==, approx(0));
    TEST(s.getVarVal("x"), ==, approx(1));

    TEST(s.getObjVal(), ==, approx(3));
  }
  {
    GurobiSolver s(shared::optim::MAX);

    size_t col1 = s.addCol("x", shared::optim::BIN, 0);
    size_t col2 = s.addCol("y", shared::optim::BIN, 0);
    size_t col3 = s.addCol("z", shared::optim::BIN, 0);

    s.update();

    s.setObjCoef(col1, 1);
    s.setObjCoef(col2, 1);
    s.setObjCoef("z", 2);

    s.update();

    TEST(s.getVarByName("x"), ==, 0);
    TEST(s.getVarByName("y"), ==, 1);
    TEST(s.getVarByName("z"), ==, 2);

    size_t row1 = s.addRow("constr1", 4, shared::optim::UP);
    s.addColToRow(row1, col1, 1);
    s.addColToRow(row1, col2, 2);
    s.addColToRow(row1, col3, 3);

    size_t row2 = s.addRow("constr2", 1, shared::optim::LO);
    s.addColToRow(row2, col1, 1);
    s.addColToRow(row2, col2, 1);

    s.update();

    TEST(s.getConstrByName("constr1"), ==, 0);
    TEST(s.getConstrByName("constr2"), ==, 1);

    auto ret = s.solve();

    TEST(ret, ==, shared::optim::OPTIM);

    TEST(s.getVarVal("x"), ==, approx(1));
    TEST(s.getVarVal("y"), ==, approx(0));
    TEST(s.getVarVal("x"), ==, approx(1));

    TEST(s.getObjVal(), ==, approx(3));
  }

#endif
}
