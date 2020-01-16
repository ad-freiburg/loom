// Copyright 2016
// Author: Patrick Brosi

#include <cassert>
#include <string>
#include "shared/tests/GurobiSolverTest.h"

#ifdef GUROBI_FOUND

#include "shared/optim/GurobiSolver.h"

using shared::optim::GurobiSolver;

#endif

// _____________________________________________________________________________
void GurobiSolverTest::run() {
#ifdef GUROBI_FOUND
  GurobiSolver s;
#endif
}
