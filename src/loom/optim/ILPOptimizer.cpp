// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <thread>

#include "loom/optim/ILPOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "shared/optim/ILPSolvProv.h"
#include "shared/rendergraph/OrderCfg.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using shared::linegraph::Line;
using shared::optim::ILPSolver;
using shared::rendergraph::HierarOrderCfg;
using util::DEBUG;
using util::ERROR;
using util::WARN;

// _____________________________________________________________________________
double ILPOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                  HierarOrderCfg* hc, size_t depth,
                                  OptResStats& stats) const {
  // avoid building the entire ILP for small search sizes
  if (solutionSpaceSize(g) < 500) {
    return _exhausOpt.optimizeComp(og, g, hc, depth + 1, stats);
  }

  LOGTO(DEBUG, std::cerr) << "Creating ILP problem... ";
  T_START(build);
  auto lp = createProblem(og, g);
  double buildT = T_STOP(build);
  LOGTO(DEBUG, std::cerr) << " .. done";

  if (lp->getNumVars() > static_cast<int>(stats.maxNumColsPerComp))
    stats.maxNumColsPerComp = lp->getNumVars();
  if (lp->getNumConstrs() > static_cast<int>(stats.maxNumRowsPerComp))
    stats.maxNumRowsPerComp = lp->getNumConstrs();

  if (_cfg->MPSOutputPath.size()) {
    lp->writeMps(_cfg->MPSOutputPath);
  }

  if (_cfg->ilpTimeLimit >= 0) lp->setTimeLim(_cfg->ilpTimeLimit);
  if (_cfg->ilpNumThreads != 0) lp->setNumThreads(_cfg->ilpNumThreads);

  LOGTO(DEBUG, std::cerr) << "Solving ILP problem...";

  T_START(solve);

  auto status = lp->solve();

  double solveT = T_STOP(solve);

  if (status == shared::optim::SolveType::INF) {
    LOG(WARN)
        << "No solution found for ILP problem (most likely because of a time "
           "limit)!";

    // apply null optimizer to guarantee a valid solution (the input solution)
    _nullOpt.optimizeComp(og, g, hc, depth + 1, stats);
  } else {
    LOGTO(DEBUG, std::cerr) << "(stats) ILP obj = " << lp->getObjVal();
    LOGTO(DEBUG, std::cerr) << "(stats) ILP build time = " << buildT << " ms";
    LOGTO(DEBUG, std::cerr) << "(stats) ILP solve time = " << solveT << " ms";
    if (status == shared::optim::SolveType::OPTIM)
      LOGTO(DEBUG, std::cerr) << "(stats) (which is optimal)";

    getConfigurationFromSolution(lp, hc, g);
  }

  delete lp;

  return solveT;
}

// _____________________________________________________________________________
int ILPOptimizer::getCrossingPenaltySameSeg(const OptNode* n) const {
  return _scorer.getCrossingPenSameSeg(n);
}

// _____________________________________________________________________________
int ILPOptimizer::getCrossingPenaltyDiffSeg(const OptNode* n) const {
  return _scorer.getCrossingPenDiffSeg(n);
}

// _____________________________________________________________________________
int ILPOptimizer::getSeparationPenalty(const OptNode* n) const {
  // double the value because we only count a separation once for each pair!
  return _scorer.getSeparationPen(n) * 1;
}

// _____________________________________________________________________________
void ILPOptimizer::getConfigurationFromSolution(
    ILPSolver* lp, HierarOrderCfg* hc, const std::set<OptNode*>& g) const {
  // build name index for faster lookup

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto lnEdgPart : e->pl().lnEdgParts) {
        if (lnEdgPart.wasCut) continue;
        for (size_t tp = 0; tp < e->pl().getCardinality(); tp++) {
          bool found = false;
          for (auto lo : e->pl().getLines()) {
            std::string varName = getILPVarName(e, lo.line, tp);

            double val = lp->getVarVal(varName);

            if (val > 0.5) {
              for (auto rel : lo.relatives) {
                // retrieve the original route pos
                size_t p = lnEdgPart.lnEdg->pl().linePos(rel);

                if (!(lnEdgPart.dir ^ e->pl().lnEdgParts.front().dir)) {
                  (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].insert(
                      (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].begin(), p);
                } else {
                  (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].push_back(p);
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
  ILPSolver* lp = shared::optim::getSolver(_cfg->ilpSolver, shared::optim::MIN);

  // for every segment s, we define |L(s)|^2 decision variables x_slp
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      // get string repr of lineedge part

      int rowA = lp->getNumConstrs();

      for (size_t p = 0; p < e->pl().getCardinality(); p++) {
        std::stringstream rowName;

        rowName << "sum(" << e->pl().getStrRepr() << ",p=" << p << ")";
        lp->addRow(rowName.str(), 1, shared::optim::FIX);
      }

      for (auto l : e->pl().getLines()) {
        // constraint: the sum of all x_slp over p must be 1 for equal sl
        std::stringstream rowName;
        rowName << "sum(" << e->pl().getStrRepr() << ",l=" << l.line << ")";

        int row = lp->addRow(rowName.str(), 1, shared::optim::FIX);

        for (size_t p = 0; p < e->pl().getCardinality(); p++) {
          std::string varName = getILPVarName(e, l.line, p);
          int curCol = lp->addCol(varName, shared::optim::BIN, 0);

          lp->addColToRow(row, curCol, 1);
          lp->addColToRow(rowA + p, curCol, 1);
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
  UNUSED(og);
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

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->pl().getStrRepr() << ","
             << segmentB->pl().getStrRepr() << "," << linepair.first.line << "("
             << linepair.first.line->id() << ")," << linepair.second.line << "("
             << linepair.second.line->id() << ")," << node << ")";

          int decisionVar = lp->addCol(
              ss.str(), shared::optim::BIN,
              getCrossingPenaltySameSeg(node)
                  // multiply the penalty with the number of collapsed lines!
                  * (linepair.first.relatives.size()) *
                  (linepair.second.relatives.size()));

          // introduce dec var for sep
          std::stringstream sss;
          sss << "x||_dec(" << segmentA->pl().getStrRepr() << ","
              << segmentB->pl().getStrRepr() << "," << linepair.first.line
              << "(" << linepair.first.line->id() << "),"
              << linepair.second.line << "(" << linepair.second.line->id()
              << ")," << node << ")";

          int decisionVarSep = 0;
          if (separationOpt()) {
            decisionVarSep = lp->addCol(sss.str(), shared::optim::BIN,
                                        getSeparationPenalty(node));
          }

          for (PosComPair poscomb :
               getPositionCombinations(segmentA, segmentB)) {
            if (crosses(node, segmentA, segmentB, poscomb)) {
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

              std::stringstream ss;
              ss << "dec_sum(" << segmentA->pl().getStrRepr() << ","
                 << segmentB->pl().getStrRepr() << "," << linepair.first.line
                 << "," << linepair.second.line << "pa=" << poscomb.first.first
                 << ",pb=" << poscomb.second.first
                 << ",pa'=" << poscomb.first.second
                 << ",pb'=" << poscomb.second.second << ",n=" << node << ")";

              int row = lp->addRow(ss.str(), 3, shared::optim::UP);

              lp->addColToRow(row, lineAinAatP, 1);
              lp->addColToRow(row, lineBinAatP, 1);
              lp->addColToRow(row, lineAinBatP, 1);
              lp->addColToRow(row, lineBinBatP, 1);
              lp->addColToRow(row, decisionVar, -1);
            }

            if (separationOpt() && separates(poscomb)) {
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

              std::stringstream ss;
              ss << "dec_sum_sep(" << segmentA->pl().getStrRepr() << ","
                 << segmentB->pl().getStrRepr() << "," << linepair.first.line
                 << "," << linepair.second.line << "pa=" << poscomb.first.first
                 << ",pb=" << poscomb.second.first
                 << ",pa'=" << poscomb.first.second
                 << ",pb'=" << poscomb.second.second << ",n=" << node << ")";

              int row = lp->addRow(ss.str(), 3, shared::optim::UP);

              lp->addColToRow(row, lineAinAatP, 1);
              lp->addColToRow(row, lineBinAatP, 1);
              lp->addColToRow(row, lineAinBatP, 1);
              lp->addColToRow(row, lineBinBatP, 1);
              lp->addColToRow(row, decisionVarSep, -1);
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
  UNUSED(og);
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

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->pl().getStrRepr() << ","
             << segments.first->pl().getStrRepr()
             << segments.second->pl().getStrRepr() << "," << linepair.first.line
             << "(" << linepair.first.line->id() << ")," << linepair.second.line
             << "(" << linepair.second.line->id() << ")," << node << ")";

          int decisionVar = lp->addCol(
              ss.str(), shared::optim::BIN,
              getCrossingPenaltyDiffSeg(node)
                  // multiply the penalty with the number of collapsed lines!
                  * (linepair.first.relatives.size()) *
                  (linepair.second.relatives.size()));

          for (PosCom poscomb : getPositionCombinations(segmentA)) {
            if (crosses(node, segmentA, segments, poscomb)) {
              int lineAinAatP = lp->getVarByName(
                  getILPVarName(segmentA, linepair.first.line, poscomb.first));
              int lineBinAatP = lp->getVarByName(getILPVarName(
                  segmentA, linepair.second.line, poscomb.second));

              assert(lineAinAatP > -1);
              assert(lineBinAatP > -1);

              std::stringstream ss;
              ss << "dec_sum(" << segmentA->pl().getStrRepr() << ","
                 << segments.first->pl().getStrRepr()
                 << segments.second->pl().getStrRepr() << ","
                 << linepair.first.line << "," << linepair.second.line
                 << "pa=" << poscomb.first << ",pb=" << poscomb.second
                 << ",n=" << node << ")";

              int row = lp->addRow(ss.str(), 1, shared::optim::UP);

              lp->addColToRow(row, lineAinAatP, 1);
              lp->addColToRow(row, lineBinAatP, 1);
              lp->addColToRow(row, decisionVar, -1);
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
  // TODO: this is already implemented in Optimizer!!
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
bool ILPOptimizer::separationOpt() const { return _scorer.optimizeSep(); }
