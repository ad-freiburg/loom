// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <thread>
#include "shared/optim/GurobiSolver.h"
#include "transitmap/graph/OrderCfg.h"
#include "transitmap/optim/ILPOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using shared::linegraph::Line;
using shared::optim::ILPSolver;
using shared::optim::GurobiSolver;

// _____________________________________________________________________________
int ILPOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                               HierarOrderCfg* hc, size_t depth) const {
  LOG(DEBUG) << "Creating ILP problem... ";
  auto lp = createProblem(og, g);
  LOG(DEBUG) << " .. done";

  // LOG(INFO) << "(stats) ILP has " << glp_get_num_cols(lp) << " cols and "
  // << glp_get_num_rows(lp) << " rows.";

  // if (!_cfg->glpkMPSOutputPath.empty()) {
  // LOG(DEBUG) << "Writing ILP as .mps to '" << _cfg->glpkMPSOutputPath << "'";
  // glp_write_mps(lp, GLP_MPS_FILE, 0, _cfg->glpkMPSOutputPath.c_str());
  // }

  LOG(DEBUG) << "Solving problem...";

  // if (!_cfg->externalSolver.empty()) {
  // preSolveCoinCbc(lp);
  // }

  std::chrono::high_resolution_clock::time_point t1 =
      std::chrono::high_resolution_clock::now();

  lp->solve();

  std::chrono::high_resolution_clock::time_point t2 =
      std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();

  if (_cfg->externalSolver.empty()) {
    LOG(INFO) << " === Solve done in " << duration << " ms ===";
  }

  LOG(INFO) << "(stats) ILP obj = " << lp->getObjVal();

  if (!_cfg->glpkSolutionOutputPath.empty()) {
    // LOG(DEBUG) << "Writing ILP full solution to '"
    // << _cfg->glpkSolutionOutputPath << "'";
    // glp_print_mip(lp, _cfg->glpkSolutionOutputPath.c_str());
  }

  getConfigurationFromSolution(lp, hc, g);

  return 0;
}

// _____________________________________________________________________________
int ILPOptimizer::getCrossingPenaltySameSeg(const OptNode* n) const {
  return _scorer->getCrossingPenaltySameSeg(n->pl().node);
}

// _____________________________________________________________________________
int ILPOptimizer::getCrossingPenaltyDiffSeg(const OptNode* n) const {
  return _scorer->getCrossingPenaltyDiffSeg(n->pl().node);
}

// _____________________________________________________________________________
int ILPOptimizer::getSplittingPenalty(const OptNode* n) const {
  // double the value because we only count a splitting once for each pair!
  return _scorer->getSplittingPenalty(n->pl().node) * 1;
}

// _____________________________________________________________________________
void ILPOptimizer::getConfigurationFromSolution(
    ILPSolver* lp, HierarOrderCfg* hc, const std::set<OptNode*>& g) const {
  // build name index for faster lookup
  // glp_create_index(lp);

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto etgp : e->pl().etgs) {
        if (etgp.wasCut) continue;
        for (size_t tp = 0; tp < e->pl().getCardinality(); tp++) {
          bool found = false;
          for (auto lo : e->pl().getLines()) {
            std::string varName = getILPVarName(e, lo.line, tp);

            // size_t i = glp_find_col(lp, varName.c_str());
            // assert(i > 0);
            // double val = glp_mip_col_val(lp, i);

            double val = lp->getVarVal(varName);

            if (val > 0.5) {
              for (auto rel : lo.relatives) {
                // retrieve the original route pos
                size_t p = etgp.etg->pl().linePos(rel);

                if (!(etgp.dir ^ e->pl().etgs.front().dir)) {
                  (*hc)[etgp.etg][etgp.order].insert(
                      (*hc)[etgp.etg][etgp.order].begin(), p);
                } else {
                  (*hc)[etgp.etg][etgp.order].push_back(p);
                }
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
ILPSolver* ILPOptimizer::createProblem(OptGraph* og,
                                       const std::set<OptNode*>& g) const {
  // glp_prob* lp = glp_create_prob();

  // glp_set_prob_name(lp, "edgeorder");
  // glp_set_obj_dir(lp, GLP_MIN);

  // VariableMatrix vm;

  auto lp = new GurobiSolver(shared::optim::MIN);

  // for every segment s, we define |L(s)|^2 decision variables x_slp
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      // get string repr of etg

      // size_t newCols = e->pl().getCardinality() * e->pl().getCardinality();
      // size_t cols = glp_add_cols(lp, newCols);
      size_t i = 0;
      // size_t rowA = glp_add_rows(lp, e->pl().getCardinality());

      int cols = lp->getNumVars();
      int rowA = lp->getNumConstrs();

      for (size_t p = 0; p < e->pl().getCardinality(); p++) {
        std::stringstream rowName;

        rowName << "sum(" << e->pl().getStrRepr() << ",p=" << p << ")";
        // glp_set_row_name(lp, rowA + p, rowName.str().c_str());
        // glp_set_row_bnds(lp, rowA + p, GLP_FX, 1, 1);
        lp->addRow(rowName.str(), 1, shared::optim::FIX);
      }

      for (auto l : e->pl().getLines()) {
        // constraint: the sum of all x_slp over p must be 1 for equal sl
        // size_t row = glp_add_rows(lp, 1);
        std::stringstream rowName;
        rowName << "sum(" << e->pl().getStrRepr() << ",l=" << l.line << ")";
        // glp_set_row_name(lp, row, rowName.str().c_str());
        // glp_set_row_bnds(lp, row, GLP_FX, 1, 1);

        int row = lp->addRow(rowName.str(), 1, shared::optim::FIX);

        for (size_t p = 0; p < e->pl().getCardinality(); p++) {
          std::string varName = getILPVarName(e, l.line, p);
          // size_t curCol = cols + i;
          // glp_set_col_name(lp, curCol, varName.c_str());
          int curCol = lp->addCol(varName, shared::optim::BIN, 0);

          // binary variable € {0,1}
          // glp_set_col_kind(lp, curCol, GLP_BV);

          // vm.addVar(row, curCol, 1);
          // vm.addVar(rowA + p, curCol, 1);

          lp->addColToRow(row, curCol, 1);
          lp->addColToRow(rowA + p, curCol, 1);

          i++;
        }
      }
    }
  }

  lp->update();

  writeSameSegConstraints(og, g, lp);
  writeDiffSegConstraints(og, g, lp);

  return lp;
}

// _____________________________________________________________________________
void ILPOptimizer::writeSameSegConstraints(OptGraph* og,
                                           const std::set<OptNode*>& g,
                                           ILPSolver* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->getAdjList()) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA)) {
        // iterate over all edges this
        // pair traverses to _TOGETHER_
        // (its possible that there are multiple edges if a line continues
        //  in more then 1 segment)
        for (OptEdge* segmentB : getEdgePartners(node, segmentA, linepair)) {
          if (processed.find(segmentB) != processed.end()) continue;
          // try all position combinations
          // size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->pl().getStrRepr() << ","
             << segmentB->pl().getStrRepr() << "," << linepair.first.line << "("
             << linepair.first.line->id() << ")," << linepair.second.line << "("
             << linepair.second.line->id() << ")," << node << ")";
          // glp_set_col_name(lp, decisionVar, ss.str().c_str());
          // glp_set_col_kind(lp, decisionVar, GLP_BV);
          // glp_set_obj_coef(
          // lp, decisionVar,
          // getCrossingPenaltySameSeg(node)
          // // multiply the penalty with the number of collapsed lines!
          // * (linepair.first.relatives.size()) *
          // (linepair.second.relatives.size()));

          int decisionVar = lp->addCol(
              ss.str(), shared::optim::BIN,
              getCrossingPenaltySameSeg(node)
                  // multiply the penalty with the number of collapsed lines!
                  * (linepair.first.relatives.size()) *
                  (linepair.second.relatives.size()));

          for (PosComPair poscomb :
               getPositionCombinations(segmentA, segmentB)) {
            if (crosses(og, node, segmentA, segmentB, poscomb)) {
              // size_t lineAinAatP =
              // glp_find_col(lp, getILPVarName(segmentA, linepair.first.line,
              // poscomb.first.first)
              // .c_str());
              // size_t lineBinAatP =
              // glp_find_col(lp, getILPVarName(segmentA, linepair.second.line,
              // poscomb.second.first)
              // .c_str());
              // size_t lineAinBatP =
              // glp_find_col(lp, getILPVarName(segmentB, linepair.first.line,
              // poscomb.first.second)
              // .c_str());
              // size_t lineBinBatP =
              // glp_find_col(lp, getILPVarName(segmentB, linepair.second.line,
              // poscomb.second.second)
              // .c_str());

              int lineAinAatP = lp->getVarByName(getILPVarName(
                  segmentA, linepair.first.line, poscomb.first.first));
              int lineBinAatP = lp->getVarByName(getILPVarName(
                  segmentA, linepair.second.line, poscomb.second.first));
              int lineAinBatP = lp->getVarByName(getILPVarName(
                  segmentB, linepair.first.line, poscomb.first.second));
              int lineBinBatP = lp->getVarByName(getILPVarName(
                  segmentB, linepair.second.line, poscomb.second.second));

              assert(lineAinAatP > -1);
              assert(lineAinBatP > -1);
              assert(lineBinAatP > -1);
              assert(lineBinBatP > -1);

              // size_t row = glp_add_rows(lp, 1);
              std::stringstream ss;
              ss << "dec_sum(" << segmentA->pl().getStrRepr() << ","
                 << segmentB->pl().getStrRepr() << "," << linepair.first.line
                 << "," << linepair.second.line << "pa=" << poscomb.first.first
                 << ",pb=" << poscomb.second.first
                 << ",pa'=" << poscomb.first.second
                 << ",pb'=" << poscomb.second.second << ",n=" << node << ")";
              // glp_set_row_name(lp, row, ss.str().c_str());
              // glp_set_row_bnds(lp, row, GLP_UP, 0, 3);

              int row = lp->addRow(ss.str(), 3, shared::optim::UP);

              lp->addColToRow(row, lineAinAatP, 1);
              lp->addColToRow(row, lineBinAatP, 1);
              lp->addColToRow(row, lineAinBatP, 1);
              lp->addColToRow(row, lineBinBatP, 1);
              lp->addColToRow(row, decisionVar, -1);

              // vm->addVar(row, lineAinAatP, 1);
              // vm->addVar(row, lineBinAatP, 1);
              // vm->addVar(row, lineAinBatP, 1);
              // vm->addVar(row, lineBinBatP, 1);
              // vm->addVar(row, decisionVar, -1);
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ILPOptimizer::writeDiffSegConstraints(OptGraph* og,
                                           const std::set<OptNode*>& g,
                                           ILPSolver* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->getAdjList()) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA)) {
        for (EdgePair segments :
             getEdgePartnerPairs(node, segmentA, linepair)) {
          // try all position combinations
          // size_t decisionVar = glp_add_cols(lp, 1);

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->pl().getStrRepr() << ","
             << segments.first->pl().getStrRepr()
             << segments.second->pl().getStrRepr() << "," << linepair.first.line
             << "(" << linepair.first.line->id() << ")," << linepair.second.line
             << "(" << linepair.second.line->id() << ")," << node << ")";
          // glp_set_col_name(lp, decisionVar, ss.str().c_str());
          // glp_set_col_kind(lp, decisionVar, GLP_BV);
          // glp_set_obj_coef(
              // lp, decisionVar,
              // getCrossingPenaltyDiffSeg(node)
                  // // multiply the penalty with the number of collapsed lines!
                  // * (linepair.first.relatives.size()) *
                  // (linepair.second.relatives.size()));

          int decisionVar = lp->addCol(ss.str(), shared::optim::BIN,
              getCrossingPenaltyDiffSeg(node)
                  // multiply the penalty with the number of collapsed lines!
                  * (linepair.first.relatives.size()) *
                  (linepair.second.relatives.size()));

          for (PosCom poscomb : getPositionCombinations(segmentA)) {
            if (crosses(og, node, segmentA, segments, poscomb)) {
              // size_t lineAinAatP = glp_find_col(
                  // lp,
                  // getILPVarName(segmentA, linepair.first.line, poscomb.first)
                      // .c_str());
              // size_t lineBinAatP = glp_find_col(
                  // lp,
                  // getILPVarName(segmentA, linepair.second.line, poscomb.second)
                      // .c_str());

              int lineAinAatP = lp->getVarByName(
                  getILPVarName(segmentA, linepair.first.line, poscomb.first)
                      );
              int lineBinAatP = lp->getVarByName(
                  getILPVarName(segmentA, linepair.second.line, poscomb.second)
                      );

              assert(lineAinAatP > -1);
              assert(lineBinAatP > -1);

              // size_t row = glp_add_rows(lp, 1);
              std::stringstream ss;
              ss << "dec_sum(" << segmentA->pl().getStrRepr() << ","
                 << segments.first->pl().getStrRepr()
                 << segments.second->pl().getStrRepr() << ","
                 << linepair.first.line << "," << linepair.second.line
                 << "pa=" << poscomb.first << ",pb=" << poscomb.second
                 << ",n=" << node << ")";
              // glp_set_row_name(lp, row, ss.str().c_str());
              // glp_set_row_bnds(lp, row, GLP_UP, 0, 1);

              int row = lp->addRow(ss.str(), 1, shared::optim::UP);

              lp->addColToRow(row, lineAinAatP, 1);
              lp->addColToRow(row, lineBinAatP, 1);
              lp->addColToRow(row, decisionVar, -1);

              // vm->addVar(row, lineAinAatP, 1);
              // vm->addVar(row, lineBinAatP, 1);
              // vm->addVar(row, decisionVar, -1);
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
std::vector<PosComPair> ILPOptimizer::getPositionCombinations(
    OptEdge* a, OptEdge* b) const {
  std::vector<PosComPair> ret;
  for (size_t posLineAinA = 0; posLineAinA < a->pl().getCardinality();
       posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < a->pl().getCardinality();
         posLineBinA++) {
      if (posLineAinA == posLineBinA) continue;

      for (size_t posLineAinB = 0; posLineAinB < b->pl().getCardinality();
           posLineAinB++) {
        for (size_t posLineBinB = 0; posLineBinB < b->pl().getCardinality();
             posLineBinB++) {
          if (posLineAinB == posLineBinB) continue;

          ret.push_back(PosComPair(PosCom(posLineAinA, posLineAinB),
                                   PosCom(posLineBinA, posLineBinB)));
        }
      }
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::vector<PosCom> ILPOptimizer::getPositionCombinations(OptEdge* a) const {
  std::vector<PosCom> ret;
  for (size_t posLineAinA = 0; posLineAinA < a->pl().getCardinality();
       posLineAinA++) {
    for (size_t posLineBinA = 0; posLineBinA < a->pl().getCardinality();
         posLineBinA++) {
      if (posLineAinA == posLineBinA) continue;
      ret.push_back(PosCom(posLineAinA, posLineBinA));
    }
  }
  return ret;
}

// _____________________________________________________________________________
std::string ILPOptimizer::getILPVarName(OptEdge* seg, const Line* r,
                                        size_t p) const {
  std::stringstream varName;
  varName << "x_(" << seg->pl().getStrRepr() << ",l=" << r << ",p=" << p << ")";
  return varName.str();
}

// _____________________________________________________________________________
// void ILPOptimizer::solveProblem(glp_prob* lp) const {
  // glp_iocp params;
  // glp_init_iocp(&params);
  // params.presolve = GLP_ON;
  // params.binarize = GLP_ON;
  // params.ps_tm_lim = _cfg->glpkPSTimeLimit;
  // params.tm_lim = _cfg->glpkTimeLimit;
  // if (_cfg->externalSolver.empty()) {
    // params.fp_heur = _cfg->useGlpkFeasibilityPump ? GLP_ON : GLP_OFF;
    // params.ps_heur = _cfg->useGlpkProximSearch ? GLP_ON : GLP_OFF;
  // }

  // glp_intopt(lp, &params);
// }

// _____________________________________________________________________________
// void ILPOptimizer::preSolveCoinCbc(glp_prob* lp) const {
  // // write temporary file
  // std::string f = std::string(std::tmpnam(0)) + ".mps";
  // std::string outf = std::string(std::tmpnam(0)) + ".sol";
  // glp_write_mps(lp, GLP_MPS_FILE, 0, f.c_str());
  // LOG(INFO) << "Calling external solver...";

  // std::chrono::high_resolution_clock::time_point t1 =
      // std::chrono::high_resolution_clock::now();

  // std::string cmd = _cfg->externalSolver;
  // util::replaceAll(cmd, "{INPUT}", f);
  // util::replaceAll(cmd, "{OUTPUT}", outf);
  // util::replaceAll(cmd, "{THREADS}",
                   // util::toString(std::thread::hardware_concurrency()));
  // LOG(INFO) << "Cmd: '" << cmd << "'";
  // int r = system(cmd.c_str());

  // std::chrono::high_resolution_clock::time_point t2 =
      // std::chrono::high_resolution_clock::now();
  // auto duration =
      // std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
  // LOG(INFO) << " === External solve done (ret=" << r << ") in " << duration
            // << " ms ===";
  // LOG(INFO) << "Parsing solution...";

  // std::ifstream fin;
  // fin.open(outf.c_str());
  // std::string line;

  // // skip first line
  // std::getline(fin, line);

  // while (std::getline(fin, line)) {
    // std::istringstream iss(line);
    // int number;
    // string name;
    // double value;

    // iss >> number;

    // // could not read line number, is missing, which is fine for us
    // if (iss.fail()) iss.clear();

    // iss >> name;
    // iss >> value;

    // int intVal = value;

    // size_t col = glp_find_col(lp, name.c_str());
    // if (col != 0) {
      // glp_set_col_bnds(lp, col, GLP_FX, intVal, intVal);
    // }
  // }
// }
