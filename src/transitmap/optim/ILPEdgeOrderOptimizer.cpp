// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include <fstream>
#include "transitmap/graph/OrderingConfiguration.h"
#include "transitmap/output/OgrOutput.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "util/geo/Geo.h"
#include "util/log/Log.h"

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
        for (size_t tp = 0; tp < etgp.etg->getCardinality(true); tp++) {
          bool found = false;
          for (size_t p = 0; p < etgp.etg->getCardinality(); p++) {
            auto r = (*etgp.etg->getTripsUnordered())[p];
            if (r.route->relativeTo()) continue;

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

  expandRelatives(c, g.getGraph());
}

// _____________________________________________________________________________
glp_prob* ILPEdgeOrderOptimizer::createProblem(const OptGraph& g) const {
  glp_prob* lp = glp_create_prob();

  glp_set_prob_name(lp, "edgeorder_impr");
  glp_set_obj_dir(lp, GLP_MIN);

  VariableMatrix vm;
  std::set<OptEdge*> processed;

  for (OptNode* n : g.getNodes()) {
    for (OptEdge* e : n->adjListOut) {
      // the first stored etg is always the ref
      graph::Edge* etg = e->etgs[0].etg;

      // constraint: the sum of all x_sl<=p over the set of lines
      // must be p+1

      size_t rowA = glp_add_rows(lp, etg->getCardinality(true));
      for (size_t p = 0; p < etg->getCardinality(true); p++) {
        std::stringstream rowName;
        rowName << "sum(" << e->getStrRepr() << ",<=" << p << ")";
        glp_set_row_name(lp, rowA + p, rowName.str().c_str());
        glp_set_row_bnds(lp, rowA + p, GLP_FX, p + 1, p + 1);
      }

      for (auto r : *etg->getTripsUnordered()) {
        if (r.route->relativeTo()) continue;
        for (size_t p = 0; p < etg->getCardinality(true); p++) {
          std::stringstream varName;
          varName << "x_(" << e->getStrRepr() << ",l=" << r.route << ",p<=" << p
                  << ")";
          size_t curCol = glp_add_cols(lp, 1);
          glp_set_col_name(lp, curCol, varName.str().c_str());
          glp_set_col_kind(lp, curCol, GLP_BV);

          // coefficients for constraint from above
          vm.addVar(rowA + p, curCol, 1);

          if (p > 0) {
            size_t row = glp_add_rows(lp, 1);

            std::stringstream rowName;
            rowName.str("");
            rowName << "sum(" << e->getStrRepr() << ",r=" << r.route <<  ",p<=" << p << ")";

            glp_set_row_name(lp, row, rowName.str().c_str());
            glp_set_row_bnds(lp, row, GLP_LO, 0, 1);

            vm.addVar(row, curCol, 1);
            vm.addVar(row, curCol - 1, -1);
          }
        }
      }
    }

  }

  glp_create_index(lp);

  writeCrossingOracle(g, &vm, lp);
  writeDiffSegConstraintsImpr(g, &vm, lp);

  int* ia = 0;
  int* ja = 0;
  double* res = 0;
  vm.getGLPKArrs(&ia, &ja, &res);

  glp_load_matrix(lp, vm.getNumVars(), ia, ja, res);

  delete[](ia);
  delete[](ja);
  delete[](res);

  return lp;
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::writeCrossingOracle(const OptGraph& g,
                                                VariableMatrix* vm,
                                                glp_prob* lp) const {
  // do everything iteratively, otherwise it would be unreadable

  size_t m = 0;

  // introduce crossing constraint variables
  for (OptNode* node : g.getNodes()) {
    for (OptEdge* segment : node->adjListOut) {
      if (segment->etgs[0].etg->getCardinality(true) > m) {
        m = segment->etgs[0].etg->getCardinality(true);
      }

      size_t rowDistanceRangeKeeper = 0;
      size_t c = segment->etgs[0].etg->getCardinality(true);
      // constraint is only needed for segments with more than 2 lines
      if (_cfg->splittingOpt && c > 2) {
        std::stringstream rowName;
        rowName << "sum_distancorRangeKeeper(e=" << segment->getStrRepr()<< ")";
        rowDistanceRangeKeeper = glp_add_rows(lp, 1);
        glp_set_row_name(lp, rowDistanceRangeKeeper, rowName.str().c_str());
        size_t max = getLinePairs(segment).size() - (2*c - 2);
        assert(max % 2 == 0);
        max = max / 2;
        glp_set_row_bnds(lp, rowDistanceRangeKeeper, GLP_UP, max, max);
      }

      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segment)) {

        // variable to check if position of line A (first) is < than
        // position of line B (second) in segment
        size_t smallerVar = glp_add_cols(lp, 1);
        std::stringstream ss;
        ss << "x_(" << segment->getStrRepr() << "," << linepair.first << "<"
           << linepair.second << ")";
        glp_set_col_name(lp, smallerVar, ss.str().c_str());
        glp_set_col_kind(lp, smallerVar, GLP_BV);
      }

      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segment, true)) {
        if (_cfg->splittingOpt && c > 2) {
          std::stringstream ss;
          // variable to check if distance between position of A and position
          // of B is > 1
          size_t dist1Var = glp_add_cols(lp, 1);
          ss << "x_(" << segment->getStrRepr() << "," << linepair.first << "<T>"
             << linepair.second << ")";
          glp_set_col_name(lp, dist1Var, ss.str().c_str());
          glp_set_col_kind(lp, dist1Var, GLP_BV);

          vm->addVar(rowDistanceRangeKeeper, dist1Var, 1);
        }
      }
    }
  }

  // write constraints for the A>B variable, both can never be 1...
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

        vm->addVar(row, smaller, 1);
        vm->addVar(row, bigger, 1);
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
        size_t rowSmallerThan = glp_add_rows(lp, 1);
        glp_set_row_name(lp, rowSmallerThan, rowName.str().c_str());
        glp_set_row_bnds(lp, rowSmallerThan, GLP_LO, 0, 0);

        std::stringstream ss;
        ss << "x_(" << segment->getStrRepr() << "," << linepair.first << "<"
           << linepair.second << ")";

        size_t decVar = glp_find_col(lp, ss.str().c_str());
        assert(decVar > 0);

        vm->addVar(rowSmallerThan, decVar, m);

        for (size_t p = 0; p < segment->etgs[0].etg->getCardinality(true); ++p) {
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

          std::stringstream rowNameSmallerThan;
          rowNameSmallerThan << "sum_crossor(e=" << segment->getStrRepr()
                  << ",A=" << linepair.first << ",B=" << linepair.second << ")";
          size_t rowSmallerThan = glp_find_row(lp, rowNameSmallerThan.str().c_str());

          vm->addVar(rowSmallerThan, first, 1);
          vm->addVar(rowSmallerThan, second, -1);

        }
      }
    }
  }

  // sum constraint for separation
  for (OptNode* node : g.getNodes()) {
    for (OptEdge* segment : node->adjListOut) {
      for (LinePair linepair : getLinePairs(segment, true)) {
        std::stringstream rowName;

        size_t rowDistance1 = 0;
        size_t rowDistance2 = 0;
        if (_cfg->splittingOpt && segment->etgs[0].etg->getCardinality(true) > 2) {
          rowName.str("");
          rowName << "sum_distancor1(e=" << segment->getStrRepr()
                  << ",A=" << linepair.first << ",B=" << linepair.second << ")";
          rowDistance1 = glp_add_rows(lp, 1);
          glp_set_row_name(lp, rowDistance1, rowName.str().c_str());
          glp_set_row_bnds(lp, rowDistance1, GLP_UP, 0, 1);

          rowName.str("");
          rowName << "sum_distancor2(e=" << segment->getStrRepr()
                  << ",A=" << linepair.first << ",B=" << linepair.second << ")";
          rowDistance2 = glp_add_rows(lp, 1);
          glp_set_row_name(lp, rowDistance2, rowName.str().c_str());
          glp_set_row_bnds(lp, rowDistance2, GLP_UP, 0, 1);

          std::stringstream sss;
          sss << "x_(" << segment->getStrRepr() << "," << linepair.first << "<T>"
             << linepair.second << ")";

          size_t decVarDistance = glp_find_col(lp, sss.str().c_str());
          assert(decVarDistance > 0);

          vm->addVar(rowDistance1, decVarDistance, -static_cast<int>(m));
          vm->addVar(rowDistance2, decVarDistance, -static_cast<int>(m));

          for (size_t p = 0; p < segment->etgs[0].etg->getCardinality(true); ++p) {
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

            vm->addVar(rowDistance1, first, 1);
            vm->addVar(rowDistance1, second, -1);

            vm->addVar(rowDistance2, first, -1);
            vm->addVar(rowDistance2, second, 1);
          }
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
             << linepair.first << "(" << linepair.first->getId() << "),"
             << linepair.second << "(" << linepair.second->getId() << ")," << node
             << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);

          glp_set_obj_coef(lp, decisionVar,
                           getCrossingPenalty(node, _cfg->crossPenMultiSameSeg));

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
                  << ",e2=" << segmentB->getStrRepr() << ",A=" << linepair.first
                  << ",B=" << linepair.second << ",n=" << node << ")";

          size_t row = glp_add_rows(lp, 1);

          glp_set_row_name(lp, row, rowName.str().c_str());
          glp_set_row_bnds(lp, row, GLP_LO, 0, 0);

          std::stringstream rowName2;
          rowName2 << "sum_dec2(e1=" << segmentA->getStrRepr()
                  << ",e2=" << segmentB->getStrRepr() << ",A=" << linepair.first
                  << ",B=" << linepair.second << ",n=" << node << ")";

          size_t row2 = glp_add_rows(lp, 1);

          glp_set_row_name(lp, row2, rowName2.str().c_str());
          glp_set_row_bnds(lp, row2, GLP_LO, 0, 0);

          bool otherWayA =
              (segmentA->from != node) ^ segmentA->etgs.front().dir;
          bool otherWayB =
              (segmentB->from != node) ^ segmentB->etgs.front().dir;

          if (otherWayA ^ otherWayB) {
          } else {
            aSmallerBinL2 = bSmallerAinL2;
          }

          vm->addVar(row, aSmallerBinL1, -1);
          vm->addVar(row, aSmallerBinL2, 1);
          vm->addVar(row, decisionVar, 1);
          vm->addVar(row2, aSmallerBinL1, 1);
          vm->addVar(row2, aSmallerBinL2, -1);
          vm->addVar(row2, decisionVar, 1);
        }
      }

      // iterate over all unique possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA, true)) {
        // iterate over all edges this
        // pair traverses to _TOGETHER_
        // (its possible that there are multiple edges if a line continues
        //  in more then 1 segment)
        for (OptEdge* segmentB : getEdgePartners(node, segmentA, linepair)) {
          if (processed.find(segmentB) != processed.end()) continue;

          // _____________________________________
          // introduce dec var for distance 1 between lines changes
          if (_cfg->splittingOpt) {
            if (segmentA->etgs[0].etg->getCardinality(true) > 2 &&
                  segmentB->etgs[0].etg->getCardinality(true) > 2) {
              // the interesting case where the line continue together from
              // segment A to segment B and the cardinality of both A and B
              // is > 2 (that is, it is possible in A or B that the two lines
              // won't be together)
              size_t decisionVarDist1Change = glp_add_cols(lp, 1);
              std::stringstream sss;
              sss << "x_decT(" << segmentA->getStrRepr() << ","
                 << segmentA->getStrRepr() << segmentB->getStrRepr() << ","
                 << linepair.first << "(" << linepair.first->getId() << "),"
                 << linepair.second << "(" << linepair.second->getId() << ")," << node
                 << ")";
              glp_set_col_name(lp, decisionVarDist1Change, sss.str().c_str());
              glp_set_col_kind(lp, decisionVarDist1Change, GLP_BV);

              glp_set_obj_coef(lp, decisionVarDist1Change,
                                getSplittingPenalty(node, _cfg->splitPenWeight));

              size_t aNearBinL1 = 0;
              size_t aNearBinL2 = 0;

              std::stringstream aNearBL1Str;
              aNearBL1Str << "x_(" << segmentA->getStrRepr() << "," << linepair.first
                      << "<T>" << linepair.second << ")";
              aNearBinL1 = glp_find_col(lp, aNearBL1Str.str().c_str());
              assert(aNearBinL1);

              std::stringstream aNearBL2Str;
              aNearBL2Str << "x_(" << segmentB->getStrRepr() << "," << linepair.first
                      << "<T>" << linepair.second << ")";
              aNearBinL2 = glp_find_col(lp, aNearBL2Str.str().c_str());
              assert(aNearBinL2);

              std::stringstream rowTName;
              rowTName << "sum_decT(e1=" << segmentA->getStrRepr()
                      << ",e2=" << segmentB->getStrRepr() << ",A=" << linepair.first
                      << ",B=" << linepair.second << ",n=" << node << ")";

              size_t check = glp_find_row(lp, rowTName.str().c_str());
              assert(check == 0);

              size_t rowT = glp_add_rows(lp, 1);

              glp_set_row_name(lp, rowT, rowTName.str().c_str());
              glp_set_row_bnds(lp, rowT, GLP_LO, 0, 0);

              std::stringstream rowTName2;
              rowTName2 << "sum_decT2(e1=" << segmentA->getStrRepr()
                      << ",e2=" << segmentB->getStrRepr() << ",A=" << linepair.first
                      << ",B=" << linepair.second << ",n=" << node << ")";

              size_t rowT2 = glp_add_rows(lp, 1);

              glp_set_row_name(lp, rowT2, rowTName2.str().c_str());
              glp_set_row_bnds(lp, rowT2, GLP_LO, 0, 0);
              vm->addVar(rowT, aNearBinL1, -1);
              vm->addVar(rowT, aNearBinL2, 1);
              vm->addVar(rowT, decisionVarDist1Change, 1);
              vm->addVar(rowT2, aNearBinL1, 1);
              vm->addVar(rowT2, aNearBinL2, -1);
              vm->addVar(rowT2, decisionVarDist1Change, 1);
            } else if ((segmentA->etgs[0].etg->getCardinality(true) == 2) ^
                  (segmentB->etgs[0].etg->getCardinality(true) == 2)) {
              // the trivial case where one of the two segments only has
              // cardinality = 2, so the lines will always be together

              OptEdge* segment = segmentA->etgs[0].etg->getCardinality(true) != 2 ? segmentA : segmentB;

              std::stringstream aNearBStr;
              aNearBStr << "x_(" << segment->getStrRepr() << "," << linepair.first
                      << "<T>" << linepair.second << ")";
              size_t aNearB = glp_find_col(lp, aNearBStr.str().c_str());
              glp_set_obj_coef(lp, aNearB,
                                getSplittingPenalty(node, _cfg->splitPenWeight));
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::writeDiffSegConstraintsImpr(const OptGraph& g,
                                                        VariableMatrix* vm,
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

          // try all position combinations
          size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->getStrRepr() << ","
             << segments.first->getStrRepr() << segments.second->getStrRepr()
             << "," << linepair.first << "(" << linepair.first->getId() << "),"
             << linepair.second << "(" << linepair.second->getId() << ")," << node
             << ")";
          glp_set_col_name(lp, decisionVar, ss.str().c_str());
          glp_set_col_kind(lp, decisionVar, GLP_BV);

          glp_set_obj_coef(lp, decisionVar,
                           getCrossingPenalty(node, _cfg->crossPenMultiDiffSeg));

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

              vm->addVar(row, testVar, 1);
              vm->addVar(row, decisionVar, -1);

              // one cross is enough...
              break;
            }
          }
        }
      }
    }
  }
}
