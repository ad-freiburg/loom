// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "./../../log/Log.h"
#include "./ILPEdgeOrderOptimizer.h"
#include <glpk.h>

using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::optimize() {
  glp_prob* lp = createProblem();

  solveProblem(lp);

  double z = glp_get_obj_val(lp);
  double x1 = glp_get_col_prim(lp, 1);
  double x2 = glp_get_col_prim(lp, 2);

  printf("z = %g; x1 = %g; x2 = %g\n", z, x1, x2);

  glp_delete_prob(lp);
  glp_free_env();
}

// _____________________________________________________________________________
glp_prob* ILPEdgeOrderOptimizer::createProblem() const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "edgeorder");
  glp_set_obj_dir(lp, GLP_MAX);  // TODO: what does this do?

  size_t c = 0;
  int ia[10000];
  int ja[10000];
  double res[10000];

  // for every segment s, we define |L(s)|^2 decision variables x_slp
  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      if (e->getEdgeTripGeoms()->size() == 0) continue;
      assert(e->getEdgeTripGeoms()->size() == 1);

      graph::EdgeTripGeom* etg = &*(e->getEdgeTripGeoms()->begin());

      // get string repr of etg

      size_t newCols = etg->getCardinality() * etg->getCardinality();
      size_t cols = glp_add_cols(lp, newCols);

      size_t i = 0;

      size_t rowA = glp_add_rows(lp, etg->getCardinality());

      for (size_t p = 0; p < etg->getCardinality(); p++) {
        std::stringstream varName;

        varName << "sum(" << etg->getStrRepr() << ", p="
          << p << ")";
        glp_set_row_name(lp, rowA + p, varName.str().c_str());
        glp_set_row_bnds(lp, rowA + p, GLP_FX, 1, 1);
      }

      for (auto r : *etg->getTripsUnordered()) {

        // constraint: the sum of all x_slp over p must be 1 for equal sl
        size_t row = glp_add_rows(lp, 1);
        c++;
        std::stringstream varName;
        varName << "sum(" << etg->getStrRepr() << ", "
          << r.route->getShortName() << ")";
        glp_set_row_name(lp, row, varName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_FX, 1, 1);

        for (size_t p = 0; p < etg->getCardinality(); p++) {
          std::stringstream varName;
          varName << "x_(" << etg->getStrRepr() << ", "
            << r.route->getShortName() << p;
          size_t curCol = cols + i;
          glp_set_col_name(lp, curCol, varName.str().c_str());

          // binary variable € {0,1}
          glp_set_col_kind(lp, curCol, GLP_BV);

          ia[c] = row;
          ja[c] = curCol;
          res[c] = 1;

          c++;

          ia[c] = rowA + p;
          ja[c] = curCol;
          res[c] = 1;

          i++;
        }
      }
    }
  }

  glp_load_matrix(lp, c, ia, ja, res);

  return lp;
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::solveProblem(glp_prob* lp) const {
  glp_simplex(lp, 0);

}

