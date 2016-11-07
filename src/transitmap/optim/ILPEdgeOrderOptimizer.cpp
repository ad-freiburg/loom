// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "./../../log/Log.h"
#include "./ILPEdgeOrderOptimizer.h"
#include "./../graph/OrderingConfiguration.h"
#include <glpk.h>

using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::optimize() {
  glp_prob* lp = createProblem();

  // write problem for debugging...
  glp_write_mps(lp, GLP_MPS_FILE, 0, "/home/patrick/ilp");

  solveProblem(lp);

  std::cout << "ILP obj = " << glp_mip_obj_val(lp) << std::endl;

  glp_print_mip(lp, "/home/patrick/ilp.sol");

  Configuration c;
  getConfigurationFromSoluation(lp, &c);
  _g->setConfig(c);

  glp_delete_prob(lp);
  glp_free_env();
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::getConfigurationFromSoluation(glp_prob* lp,
    Configuration* c) const {
  // build name index for faster lookup
  glp_create_index(lp);

  for (graph::Node* n : *_g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      if (e->getEdgeTripGeoms()->size() == 0) continue;
      assert(e->getEdgeTripGeoms()->size() == 1);

      graph::EdgeTripGeom* etg = &*(e->getEdgeTripGeoms()->begin());

      for (size_t tp = 0; tp < etg->getCardinality(); tp++) {
        bool found = false;
        for (size_t p = 0; p < etg->getCardinality(); p++) {
          auto r = (*etg->getTripsUnordered())[p];
          std::stringstream varName;
          varName << "x_(" << etg->getStrRepr() << ",l="
            << r.route << ",p=" << tp << ")";

          size_t i = glp_find_col(lp, varName.str().c_str());
          double val = glp_mip_col_val(lp, i);

          if (val > 0.5) {
            (*c)[etg].push_back(p);
            assert(!found);  // should be assured by ILP constraints
            found = true;
          }
        }
        assert(found);
      }
    }
  }
}

// _____________________________________________________________________________
glp_prob* ILPEdgeOrderOptimizer::createProblem() const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "edgeorder");
  glp_set_obj_dir(lp, GLP_MIN);

  size_t c = 0;

  // TODO: array sizes
  int* ia = new int[1000000];
  int* ja = new int[1000000];
  double* res = new double[1000000];

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

        varName << "sum(" << etg->getStrRepr() << ",p="
          << p << ")";
        glp_set_row_name(lp, rowA + p, varName.str().c_str());
        glp_set_row_bnds(lp, rowA + p, GLP_FX, 1, 1);
      }

      for (auto r : *etg->getTripsUnordered()) {

        // constraint: the sum of all x_slp over p must be 1 for equal sl
        size_t row = glp_add_rows(lp, 1);
        std::stringstream varName;
        varName << "sum(" << etg->getStrRepr() << ",l="
          << r.route << ")";
        glp_set_row_name(lp, row, varName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_FX, 1, 1);

        for (size_t p = 0; p < etg->getCardinality(); p++) {
          std::stringstream varName;
          varName << "x_(" << etg->getStrRepr() << ",l="
            << r.route << ",p=" << p << ")";
          size_t curCol = cols + i;
          glp_set_col_name(lp, curCol, varName.str().c_str());

          // binary variable € {0,1}
          glp_set_col_kind(lp, curCol, GLP_BV);
          c++;

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

  glp_create_index(lp);

  // go into nodes and build crossing constraints for adjacent nodefronts
  for (graph::Node* n : *_g->getNodes()) {
    // if (n->getStops().size() == 0 || (*(n->getStops().begin()))->getName() != "Eschholzstraße") continue;
    for (auto& nf : n->getMainDirs()) {
      const graph::EdgeTripGeom* seg = nf.refEtg;

      std::set<gtfs::Route*> processed;

      // check all line pairs in this seg
      for (auto& toA : seg->getTripsUnordered()) {
        processed.insert(toA.route);
        for (auto& toB : seg->getTripsUnordered()) {
          if (processed.find(toB.route) != processed.end()) continue;
          if (toA.route == toB.route) continue; // dont check against itself

          std::vector<Partner> partnersA = n->getPartner(&nf, toA.route);
          std::vector<Partner> partnersB = n->getPartner(&nf, toB.route);

          const NodeFront* nfAdj = 0;

          // check for adjacent partner in A
          for (size_t i = 0; i < partnersA.size(); ++i) {
            if (n->isAdjacent(&nf, partnersA[i].front)) {
                nfAdj = partnersA[i].front;
            }
          }

          if (!nfAdj) {
            std::cout << "No adj. partner found..." << std::endl;
            continue; // no adjacent partner found
          }

          // now, check if nfa is present in partnersB also
          bool found = false;
          for (size_t i = 0; i < partnersB.size(); ++i) {
            if (partnersB[i].front == nfAdj) found = true;
          }

          if (!found) {
            std::cout << "No partner in same nf found..." << std::endl;
            continue;
          }


          const graph::EdgeTripGeom* segB = nfAdj->refEtg;

          std::cout << "-> " << toA.route << " (" << toA.route->getShortName()
            << ") and " << toB.route << " (" << toB.route->getShortName()
            << ") continue together from segment " << seg << " to segment "
            << nfAdj->refEtg << ", adjacent in node " << n
            << " (" << (n->getStops().size() > 0 ? (*(n->getStops().begin()))->getName() : "<tn>")
            << ")" << std::endl;

          size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << seg->getStrRepr() << "," << segB->getStrRepr()
            << "," << toA.route << "(" << toA.route->getShortName() << "),"
            << toB.route << "(" << toB.route->getShortName() << ")," << n << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);
          glp_set_obj_coef(lp, decisionVar, 1);

          for (size_t posLineAinA = 0; posLineAinA < seg->getCardinality(); posLineAinA++) {
            for (size_t posLineBinA = 0; posLineBinA < seg->getCardinality(); posLineBinA++) {

              if (posLineAinA == posLineBinA) continue; // already covered by constraint above

              for (size_t posLineAinB = 0; posLineAinB < segB->getCardinality(); posLineAinB++) {
                for (size_t posLineBinB = 0; posLineBinB < segB->getCardinality(); posLineBinB++) {

                  if (posLineAinB == posLineBinB) continue; // already covered by constraint above

                  if (n->crosses(*seg, *segB, posLineAinA, posLineAinB, posLineBinA, posLineBinB)) {
                    std::cout << "  crosses with pa=" << posLineAinA << ", pb=" << posLineBinA
                      << ", pa'=" << posLineAinB << ", pb'=" << posLineBinB << std::endl;
                    // get variables

                    size_t lineAinAatP = glp_find_col(lp, getILPVarName(seg, toA.route, posLineAinA).c_str());
                    size_t lineBinAatP = glp_find_col(lp, getILPVarName(seg, toB.route, posLineBinA).c_str());
                    size_t lineAinBatP = glp_find_col(lp, getILPVarName(segB, toA.route, posLineAinB).c_str());
                    size_t lineBinBatP = glp_find_col(lp, getILPVarName(segB, toB.route, posLineBinB).c_str());

                    assert(lineAinAatP > 0);
                    assert(lineAinBatP > 0);
                    assert(lineBinAatP > 0);
                    assert(lineBinBatP > 0);

                    size_t row = glp_add_rows(lp, 1);
                    std::stringstream ss;
                    ss << "dec_sum(" << seg->getStrRepr() << "," << segB->getStrRepr()
                      << "," << toA.route << "," << toB.route
                      << "pa=" << posLineAinA << ",pb=" << posLineBinA
                     << ",pa'=" << posLineAinB << ",pb'=" << posLineBinB
                      << ",n=" << n << ")";
                    glp_set_row_name(lp, row, ss.str().c_str());
                    glp_set_row_bnds(lp, row, GLP_UP, 0, 3);

                    c++;
                    ia[c] = row;
                    ja[c] = lineAinAatP;
                    res[c] = 1;

                    c++;
                    ia[c] = row;
                    ja[c] = lineBinAatP;
                    res[c] = 1;

                    c++;
                    ia[c] = row;
                    ja[c] = lineAinBatP;
                    res[c] = 1;

                    c++;
                    ia[c] = row;
                    ja[c] = lineBinBatP;
                    res[c] = 1;

                    c++;
                    ia[c] = row;
                    ja[c] = decisionVar;
                    res[c] = -1;
                  }
                }
              }

            }
          }

        }
      }
    }
  }


  glp_load_matrix(lp, c, ia, ja, res);

  delete[](ia);
  delete[](ja);
  delete[](res);

  return lp;
}

// _____________________________________________________________________________
std::string ILPEdgeOrderOptimizer::getILPVarName(const EdgeTripGeom* seg,
    const gtfs::Route* r, size_t p) const {
  std::stringstream varName;
  varName << "x_(" << seg->getStrRepr() << ",l="
    << r << ",p=" << p << ")";
  return varName.str();
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::solveProblem(glp_prob* lp) const {
  glp_iocp params;
  glp_init_iocp(&params);
  params.presolve = GLP_ON;
  params.tm_lim = 6000000;
  glp_intopt(lp, &params);
}

