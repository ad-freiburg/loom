// Copyright 2016
// Author: Patrick Brosi

#include <cassert>
#include <string>
#include <vector>
#include "shared/optim/ILPSolver.h"
#include "shared/tests/ILPSolverTest.h"
#include "util/Misc.h"
#include "util/Test.h"

using shared::optim::ILPSolver;
using util::approx;

#ifdef GUROBI_FOUND

#include "shared/optim/GurobiSolver.h"

using shared::optim::GurobiSolver;

#endif

#ifdef GLPK_FOUND

#include "shared/optim/GLPKSolver.h"

using shared::optim::GLPKSolver;

#endif

#ifdef COIN_FOUND

#include "shared/optim/COINSolver.h"

using shared::optim::COINSolver;

#endif

// _____________________________________________________________________________
void ILPSolverTest::run() {
  {
    std::vector<ILPSolver*> solvers;

#ifdef GUROBI_FOUND
    try {
      solvers.push_back(new GurobiSolver(shared::optim::MAX));
    } catch (const std::exception& e) {
    }
#endif

#ifdef GLPK_FOUND
    solvers.push_back(new GLPKSolver(shared::optim::MAX));
#endif

#ifdef COIN_FOUND
    auto s = new COINSolver(shared::optim::MAX);
    solvers.push_back(s);
#endif


    for (auto s : solvers) {
      size_t col1 = s->addCol("x", shared::optim::BIN, 1);
      size_t col2 = s->addCol("y", shared::optim::BIN, 1);
      size_t col3 = s->addCol("z", shared::optim::BIN, 2);

      s->update();

      TEST(s->getVarByName("x"), ==, 0);
      TEST(s->getVarByName("y"), ==, 1);
      TEST(s->getVarByName("z"), ==, 2);

      size_t row1 = s->addRow("constr1", 4, shared::optim::UP);
      s->addColToRow(row1, col1, 1);
      s->addColToRow(row1, col2, 2);
      s->addColToRow(row1, col3, 3);

      size_t row2 = s->addRow("constr2", 1, shared::optim::LO);
      s->addColToRow(row2, col1, 1);
      s->addColToRow(row2, col2, 1);

      s->update();

      TEST(s->getConstrByName("constr1"), ==, 0);
      TEST(s->getConstrByName("constr2"), ==, 1);

      auto ret = s->solve();

      TEST(ret, ==, shared::optim::OPTIM);

      TEST(s->getVarVal("x"), ==, approx(1));
      TEST(s->getVarVal("y"), ==, approx(0));
      TEST(s->getVarVal("z"), ==, approx(1));

      TEST(s->getObjVal(), ==, approx(3));
    }
  }
  {
    std::vector<ILPSolver*> solvers;
#ifdef GUROBI_FOUND
    try {
      solvers.push_back(new GurobiSolver(shared::optim::MAX));
    } catch (const std::exception& e) {
    }
#endif

#ifdef GLPK_FOUND
    solvers.push_back(new GLPKSolver(shared::optim::MAX));
#endif

    for (auto s : solvers) {
      size_t col1 = s->addCol("x", shared::optim::BIN, 0);
      size_t col2 = s->addCol("y", shared::optim::BIN, 0);
      size_t col3 = s->addCol("z", shared::optim::BIN, 0);

      s->update();

      s->setObjCoef(col1, 1);
      s->setObjCoef(col2, 1);
      s->setObjCoef("z", 2);

      s->update();

      TEST(s->getVarByName("x"), ==, 0);
      TEST(s->getVarByName("y"), ==, 1);
      TEST(s->getVarByName("z"), ==, 2);

      size_t row1 = s->addRow("constr1", 4, shared::optim::UP);
      s->addColToRow(row1, col1, 1);
      s->addColToRow(row1, col2, 2);
      s->addColToRow(row1, col3, 3);

      size_t row2 = s->addRow("constr2", 1, shared::optim::LO);
      s->addColToRow(row2, col1, 1);
      s->addColToRow(row2, col2, 1);

      s->update();

      TEST(s->getConstrByName("constr1"), ==, 0);
      TEST(s->getConstrByName("constr2"), ==, 1);

      auto ret = s->solve();

      TEST(ret, ==, shared::optim::OPTIM);

      TEST(s->getVarVal("x"), ==, approx(1));
      TEST(s->getVarVal("y"), ==, approx(0));
      TEST(s->getVarVal("z"), ==, approx(1));

      TEST(s->getObjVal(), ==, approx(3));
    }
  }
}
