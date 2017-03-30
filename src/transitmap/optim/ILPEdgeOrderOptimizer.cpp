// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include <fstream>
#include "./../graph/OrderingConfiguration.h"
#include "./../output/OgrOutput.h"
#include "./ILPEdgeOrderOptimizer.h"
#include "./OptGraph.h"
#include "pbutil/geo/Geo.h"
#include "pbutil/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace graph;

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::getConfigurationFromSolution(
    glp_prob* lp, Configuration* c, const OptGraph& g) const {
  // build name index for faster lookup
  glp_create_index(lp);

  for (OptNode* n : g.getNodes()) {
    for (OptEdge* e : n->adjListOut) {
      for (auto etgp : e->etgs) {
        for (size_t tp = 0; tp < etgp.etg->getCardinality(); tp++) {
          bool found = false;
          for (size_t p = 0; p < etgp.etg->getCardinality(); p++) {
            auto r = (*etgp.etg->getTripsUnordered())[p];

            // check if this route (r) switch from 0 to 1 at tp-1 and tp
            double valPrev = 0;
            std::stringstream varName;

            if (tp > 0) {
              varName << "x_(" << e->getStrRepr() << ",l=" << r.route
                      << ",p<=" << tp - 1 << ")";

              size_t i = glp_find_col(lp, varName.str().c_str());
              assert(i > 0);
              valPrev = glp_mip_col_val(lp, i);
            }

            varName.str("");
            varName << "x_(" << e->getStrRepr() << ",l=" << r.route
                    << ",p<=" << tp << ")";

            size_t i = glp_find_col(lp, varName.str().c_str());
            assert(i > 0);
            double val = glp_mip_col_val(lp, i);

            if (valPrev < 0.5 && val > 0.5) {
              // first time p is eq/greater, so it is this p
              if (!(etgp.dir ^ e->etgs.front().dir)) {
                (*c)[etgp.etg].insert((*c)[etgp.etg].begin(), p);
              } else {
                (*c)[etgp.etg].push_back(p);
              }
              assert(!found);  // should be assured by ILP constraints
              found = true;
            }
          }
          assert(found);
        }
      }
    }
  }
}

// _____________________________________________________________________________
glp_prob* ILPEdgeOrderOptimizer::createProblem(const OptGraph& g) const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "edgeorder_impr");
  glp_set_obj_dir(lp, GLP_MIN);

  size_t c = 0;

  // TODO: array sizes
  int* ia = new int[1000000];
  int* ja = new int[1000000];
  double* res = new double[1000000];

  for (OptNode* n : g.getNodes()) {
    for (OptEdge* e : n->adjListOut) {
      // the first stored etg is always the ref
      graph::Edge* etg = e->etgs[0].etg;

      // constraint: the sum of all x_sl<=p over the set of lines
      // must be p+1

      size_t rowA = glp_add_rows(lp, etg->getCardinality());
      for (size_t p = 0; p < etg->getCardinality(); p++) {
        std::stringstream varName;
        varName << "sum(" << e->getStrRepr() << ",<=" << p << ")";
        glp_set_row_name(lp, rowA + p, varName.str().c_str());
        glp_set_row_bnds(lp, rowA + p, GLP_FX, p + 1, p + 1);
      }

      for (auto r : *etg->getTripsUnordered()) {
        for (size_t p = 0; p < etg->getCardinality(); p++) {
          std::stringstream varName;
          varName << "x_(" << e->getStrRepr() << ",l=" << r.route << ",p<=" << p
                  << ")";
          size_t curCol = glp_add_cols(lp, 1);
          glp_set_col_name(lp, curCol, varName.str().c_str());
          glp_set_col_kind(lp, curCol, GLP_BV);

          // coefficients for constraint from above
          c++;

          ia[c] = rowA + p;
          ja[c] = curCol;
          res[c] = 1;

          if (p > 0) {
            size_t row = glp_add_rows(lp, 1);

            std::stringstream rowName;
            rowName.str("");
            rowName << "sum(" << e->getStrRepr() << ",p<=" << p << ")";

            glp_set_row_name(lp, row, rowName.str().c_str());
            glp_set_row_bnds(lp, row, GLP_LO, 0, 1);

            c++;

            ia[c] = row;
            ja[c] = curCol;
            res[c] = 1;

            c++;

            ia[c] = row;
            ja[c] = curCol - 1;
            res[c] = -1;
          }
        }
      }
    }
  }

  glp_create_index(lp);

  writeCrossingOracle(g, ia, ja, res, &c, lp);
  writeDiffSegConstraintsImpr(g, ia, ja, res, &c, lp);

  glp_load_matrix(lp, c, ia, ja, res);

  delete[](ia);
  delete[](ja);
  delete[](res);

  return lp;
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::writeCrossingOracle(const OptGraph& g, int* ia,
                                                int* ja, double* res, size_t* c,
                                                glp_prob* lp) const {
  // do everything iteratively, otherwise it would be unreadable

  size_t m = 0;

  // introduce crossing constraint variables
  for (OptNode* node : g.getNodes()) {
    for (OptEdge* segment : node->adjListOut) {
      if (segment->etgs[0].etg->getCardinality() > m) {
        m = segment->etgs[0].etg->getCardinality();
      }

      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segment)) {
        size_t smallerVar = glp_add_cols(lp, 1);
        std::stringstream ss;
        ss << "x_(" << segment->getStrRepr() << "," << linepair.first << "<"
           << linepair.second << ")";
        glp_set_col_name(lp, smallerVar, ss.str().c_str());
        glp_set_col_kind(lp, smallerVar, GLP_BV);
      }
    }
  }

  // write constraints for the crossing variable, both can never be 1...
  for (OptNode* node : g.getNodes()) {
    for (OptEdge* segment : node->adjListOut) {
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segment)) {
        std::stringstream ss;
        ss << "x_(" << segment->getStrRepr() << "," << linepair.first << "<"
           << linepair.second << ")";
        size_t smaller = glp_find_col(lp, ss.str().c_str());
        assert(smaller > 0);

        std::stringstream ss2;
        ss2 << "x_(" << segment->getStrRepr() << "," << linepair.second << "<"
            << linepair.first << ")";
        size_t bigger = glp_find_col(lp, ss2.str().c_str());
        assert(bigger > 0);

        std::stringstream rowName;
        rowName << "sum(" << ss.str() << "," << ss2.str() << ")";

        size_t row = glp_add_rows(lp, 1);

        glp_set_row_name(lp, row, rowName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_FX, 1, 1);

        (*c)++;
        ia[*c] = row;
        ja[*c] = smaller;
        res[*c] = 1;

        (*c)++;
        ia[*c] = row;
        ja[*c] = bigger;
        res[*c] = 1;
      }
    }
  }

  // sum constraint
  for (OptNode* node : g.getNodes()) {
    for (OptEdge* segment : node->adjListOut) {
      for (LinePair linepair : getLinePairs(segment)) {
        std::stringstream rowName;
        rowName << "sum_crossor(e=" << segment->getStrRepr()
                << ",A=" << linepair.first << ",B=" << linepair.second << ")";

        size_t row = glp_add_rows(lp, 1);

        glp_set_row_name(lp, row, rowName.str().c_str());
        glp_set_row_bnds(lp, row, GLP_LO, 0, 0);

        std::stringstream ss;
        ss << "x_(" << segment->getStrRepr() << "," << linepair.first << "<"
           << linepair.second << ")";

        size_t decVar = glp_find_col(lp, ss.str().c_str());
        assert(decVar > 0);

        (*c)++;
        ia[*c] = row;
        ja[*c] = decVar;
        res[*c] = m;

        for (size_t p = 0; p < segment->etgs[0].etg->getCardinality(); ++p) {
          std::stringstream ss;
          ss << "x_(" << segment->getStrRepr() << ",l=" << linepair.first
             << ",p<=" << p << ")";
          size_t first = glp_find_col(lp, ss.str().c_str());
          assert(first > 0);

          std::stringstream ss2;
          ss2 << "x_(" << segment->getStrRepr() << ",l=" << linepair.second
              << ",p<=" << p << ")";
          size_t second = glp_find_col(lp, ss2.str().c_str());
          assert(second > 0);

          (*c)++;
          ia[*c] = row;
          ja[*c] = first;
          res[*c] = 1;

          (*c)++;
          ia[*c] = row;
          ja[*c] = second;
          res[*c] = -1;
        }
      }
    }
  }

  // crossing constraints
  for (OptNode* node : g.getNodes()) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->adjList) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA)) {
        // iterate over all edges this
        // pair traverses to _TOGETHER_
        // (its possible that there are multiple edges if a line continues
        //  in more then 1 segment)
        for (OptEdge* segmentB : getEdgePartners(node, segmentA, linepair)) {
          if (processed.find(segmentB) != processed.end()) continue;

          // introduce dec var

          size_t decisionVar = glp_add_cols(lp, 1);
          std::stringstream ss;
          ss << "x_dec(" << segmentA->getStrRepr() << ","
             << segmentA->getStrRepr() << segmentB->getStrRepr() << ","
             << linepair.first << "(" << linepair.first->id << "),"
             << linepair.second << "(" << linepair.second->id << ")," << node
             << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);

          int coef = 2;
          coef *= node->adjList.size();

          if (node->node->getStops().size() > 0) {
            glp_set_obj_coef(lp, decisionVar, 3 * coef);
          } else {
            glp_set_obj_coef(lp, decisionVar, coef);
          }

          size_t aSmallerBinL1 = 0;
          size_t aSmallerBinL2 = 0;
          size_t bSmallerAinL2 = 0;

          std::stringstream aSmBStr;
          aSmBStr << "x_(" << segmentA->getStrRepr() << "," << linepair.first
                  << "<" << linepair.second << ")";
          aSmallerBinL1 = glp_find_col(lp, aSmBStr.str().c_str());

          std::stringstream aBgBStr;
          aBgBStr << "x_(" << segmentB->getStrRepr() << "," << linepair.first
                  << "<" << linepair.second << ")";
          aSmallerBinL2 = glp_find_col(lp, aBgBStr.str().c_str());

          std::stringstream bBgAStr;
          bBgAStr << "x_(" << segmentB->getStrRepr() << "," << linepair.second
                  << "<" << linepair.first << ")";
          bSmallerAinL2 = glp_find_col(lp, bBgAStr.str().c_str());

          assert(aSmallerBinL1 && aSmallerBinL2 && bSmallerAinL2);

          std::stringstream rowName;
          rowName << "sum_dec(e1=" << segmentA->getStrRepr()
                  << "e2=" << segmentA->getStrRepr() << ",A=" << linepair.first
                  << ",B=" << linepair.second << ")";

          size_t row = glp_add_rows(lp, 1);

          glp_set_row_name(lp, row, rowName.str().c_str());
          glp_set_row_bnds(lp, row, GLP_UP, 0, 0);

          std::stringstream rowName2;
          rowName << "sum_dec2(e1=" << segmentA->getStrRepr()
                  << "e2=" << segmentA->getStrRepr() << ",A=" << linepair.first
                  << ",B=" << linepair.second << ")";

          size_t row2 = glp_add_rows(lp, 1);

          glp_set_row_name(lp, row2, rowName2.str().c_str());
          glp_set_row_bnds(lp, row2, GLP_UP, 0, 0);

          bool otherWayA =
              (segmentA->from != node) ^ segmentA->etgs.front().dir;
          bool otherWayB =
              (segmentB->from != node) ^ segmentB->etgs.front().dir;

          if (otherWayA ^ otherWayB) {
          } else {
            aSmallerBinL2 = bSmallerAinL2;
          }

          (*c)++;
          ia[*c] = row;
          ja[*c] = aSmallerBinL1;
          res[*c] = 1;

          (*c)++;
          ia[*c] = row;
          ja[*c] = aSmallerBinL2;
          res[*c] = -1;

          (*c)++;
          ia[*c] = row;
          ja[*c] = decisionVar;
          res[*c] = -1;

          (*c)++;
          ia[*c] = row2;
          ja[*c] = aSmallerBinL1;
          res[*c] = -1;

          (*c)++;
          ia[*c] = row2;
          ja[*c] = aSmallerBinL2;
          res[*c] = 1;

          (*c)++;
          ia[*c] = row2;
          ja[*c] = decisionVar;
          res[*c] = -1;
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::writeDiffSegConstraintsImpr(const OptGraph& g,
                                                        int* ia, int* ja,
                                                        double* res, size_t* c,
                                                        glp_prob* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g.getNodes()) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->adjList) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA)) {
        for (EdgePair segments :
             getEdgePartnerPairs(node, segmentA, linepair)) {
          // if (processed.find(segmentB) != processed.end()) continue;

          // try all position combinations
          size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->getStrRepr() << ","
             << segments.first->getStrRepr() << segments.second->getStrRepr()
             << "," << linepair.first << "(" << linepair.first->id << "),"
             << linepair.second << "(" << linepair.second->id << ")," << node
             << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);

          int coef = 1;
          coef *= node->adjList.size();

          if (node->node->getStops().size() > 0) {
            glp_set_obj_coef(lp, decisionVar, 3 * coef);
          } else {
            glp_set_obj_coef(lp, decisionVar, 1 * coef);
          }

          for (PosCom poscomb : getPositionCombinations(segmentA)) {
            if (crosses(node, segmentA, segments, poscomb)) {
              size_t testVar = 0;

              if (poscomb.first > poscomb.second) {
                std::stringstream bBgAStr;
                bBgAStr << "x_(" << segmentA->getStrRepr() << ","
                        << linepair.first << "<" << linepair.second << ")";
                testVar = glp_find_col(lp, bBgAStr.str().c_str());
              } else {
                std::stringstream bBgAStr;
                bBgAStr << "x_(" << segmentA->getStrRepr() << ","
                        << linepair.second << "<" << linepair.first << ")";
                testVar = glp_find_col(lp, bBgAStr.str().c_str());
              }

              assert(testVar);

              size_t row = glp_add_rows(lp, 1);
              std::stringstream ss;
              ss << "dec_sum(" << segmentA->getStrRepr() << ","
                 << segments.first->getStrRepr()
                 << segments.second->getStrRepr() << "," << linepair.first
                 << "," << linepair.second << "pa=" << poscomb.first
                 << ",pb=" << poscomb.second << ",n=" << node << ")";
              glp_set_row_name(lp, row, ss.str().c_str());
              glp_set_row_bnds(lp, row, GLP_FX, 0, 0);

              (*c)++;
              ia[*c] = row;
              ja[*c] = testVar;
              res[*c] = 1;

              (*c)++;
              ia[*c] = row;
              ja[*c] = decisionVar;
              res[*c] = -1;

              // one cross is enough...
              break;
            }
          }
        }
      }
    }
  }
}
