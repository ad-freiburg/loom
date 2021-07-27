// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "shared/optim/ILPSolvProv.h"
#include "shared/rendergraph/OrderCfg.h"
#include "util/geo/Geo.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using shared::optim::ILPSolver;
using shared::rendergraph::HierarOrderCfg;

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::getConfigurationFromSolution(
    ILPSolver* lp, HierarOrderCfg* hc, const std::set<OptNode*>& g) const {
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      for (auto etgp : e->pl().etgs) {
        if (etgp.wasCut) continue;
        for (size_t tp = 0; tp < e->pl().getCardinality(); tp++) {
          bool found = false;

          for (auto ro : e->pl().getLines()) {
            // check if this route (r) switches from 0 to 1 at tp-1 and tp
            double valPrev = 0;
            std::stringstream varName;

            if (tp > 0) {
              varName << "x_(" << e->pl().getStrRepr() << ",l=" << ro.line
                      << ",p<=" << tp - 1 << ")";

              valPrev = lp->getVarVal(varName.str());
            }

            varName.str("");
            varName << "x_(" << e->pl().getStrRepr() << ",l=" << ro.line
                    << ",p<=" << tp << ")";

            double val = lp->getVarVal(varName.str());

            if (valPrev < 0.5 && val > 0.5) {
              // first time p is eq/greater, so it is this p
              // TODO: the latter dir is checking the 'main' direction here,
              // put this into a method in the pl()! (there, the [0] etg is
              // already taken as a ref). THIS IS A POTENTIAL BUG HERE

              for (auto rel : ro.relatives) {
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
ILPSolver* ILPEdgeOrderOptimizer::createProblem(
    OptGraph* og, const std::set<OptNode*>& g) const {
  ILPSolver* lp = shared::optim::getSolver(_cfg->ilpSolver, shared::optim::MIN);

  std::set<OptEdge*> processed;

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      // constraint: the sum of all x_sl<=p over the set of lines
      // must be p+1

      size_t rowA = lp->getNumConstrs();
      for (size_t p = 0; p < e->pl().getCardinality(); p++) {
        std::stringstream rowName;
        rowName << "sum(" << e->pl().getStrRepr() << ",<=" << p << ")";
        lp->addRow(rowName.str(), p + 1, shared::optim::FIX);
      }

      for (auto r : e->pl().getLines()) {
        for (size_t p = 0; p < e->pl().getCardinality(); p++) {
          std::stringstream varName;
          varName << "x_(" << e->pl().getStrRepr() << ",l=" << r.line
                  << ",p<=" << p << ")";
          int curCol = lp->addCol(varName.str(), shared::optim::BIN, 0);

          // coefficients for constraint from above
          lp->addColToRow(rowA + p, curCol, 1);

          if (p > 0) {
            std::stringstream rowName;
            rowName << "sum(" << e->pl().getStrRepr() << ",r=" << r.line
                    << ",p<=" << p << ")";

            int row = lp->addRow(rowName.str(), 0, shared::optim::LO);

            lp->addColToRow(row, curCol, 1);
            lp->addColToRow(row, curCol - 1, -1);
          }
        }
      }
    }
  }

  lp->update();

  writeCrossingOracle(g, lp);
  writeDiffSegConstraintsImpr(g, lp);

  return lp;
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::writeCrossingOracle(const std::set<OptNode*>& g,
                                                ILPSolver* lp) const {
  // do everything iteratively, otherwise it would be unreadable

  size_t m = 0;

  // introduce crossing constraint variables
  for (OptNode* node : g) {
    for (OptEdge* segment : node->getAdjList()) {
      if (segment->getFrom() != node) continue;
      if (segment->pl().getCardinality() > m) {
        m = segment->pl().getCardinality();
      }

      size_t rowDistanceRangeKeeper = 0;
      size_t c = segment->pl().getCardinality();
      // constraint is only needed for segments with more than 2 lines
      if (_cfg->splittingOpt && c > 2) {
        std::stringstream rowName;
        rowName << "sum_distancorRangeKeeper(e=" << segment->pl().getStrRepr()
                << ")";

        size_t max = getLinePairs(segment).size() - (2 * c - 2);
        assert(max % 2 == 0);
        max = max / 2;

        rowDistanceRangeKeeper =
            lp->addRow(rowName.str(), max, shared::optim::UP);
      }

      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segment)) {
        // variable to check if position of line A (first) is < than
        // position of line B (second) in segment
        std::stringstream ss;
        ss << "x_(" << segment->pl().getStrRepr() << "," << linepair.first.line
           << "<" << linepair.second.line << ")";

        lp->addCol(ss.str(), shared::optim::BIN, 0);
      }

      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segment, true)) {
        if (_cfg->splittingOpt && c > 2) {
          std::stringstream ss;
          // variable to check if distance between position of A and position
          // of B is > 1
          ss << "x_(" << segment->pl().getStrRepr() << ","
             << linepair.first.line << "<T>" << linepair.second.line << ")";

          size_t dist1Var = lp->addCol(ss.str(), shared::optim::BIN, 0);
          lp->addColToRow(rowDistanceRangeKeeper, dist1Var, 1);
        }
      }
    }
  }

  lp->update();

  // write constraints for the A>B variable, both can never be 1...
  for (OptNode* node : g) {
    for (OptEdge* segment : node->getAdjList()) {
      if (segment->getFrom() != node) continue;
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segment)) {
        std::stringstream ss;
        ss << "x_(" << segment->pl().getStrRepr() << "," << linepair.first.line
           << "<" << linepair.second.line << ")";

        int smaller = lp->getVarByName(ss.str());
        assert(smaller > -1);

        std::stringstream ss2;
        ss2 << "x_(" << segment->pl().getStrRepr() << ","
            << linepair.second.line << "<" << linepair.first.line << ")";

        int bigger = lp->getVarByName(ss2.str());
        assert(bigger > -1);

        std::stringstream rowName;
        rowName << "sum(" << ss.str() << "," << ss2.str() << ")";

        int row = lp->addRow(rowName.str(), 1, shared::optim::FIX);

        lp->addColToRow(row, smaller, 1);
        lp->addColToRow(row, bigger, 1);
      }
    }
  }

  // sum constraint
  for (OptNode* node : g) {
    for (OptEdge* segment : node->getAdjList()) {
      if (segment->getFrom() != node) continue;
      for (LinePair linepair : getLinePairs(segment)) {
        std::stringstream rowName;
        rowName << "sum_crossor(e=" << segment->pl().getStrRepr()
                << ",A=" << linepair.first.line << ",B=" << linepair.second.line
                << ")";
        int rowSmallerThan = lp->addRow(rowName.str(), 0, shared::optim::LO);

        std::stringstream ss;
        ss << "x_(" << segment->pl().getStrRepr() << "," << linepair.first.line
           << "<" << linepair.second.line << ")";

        int decVar = lp->getVarByName(ss.str());
        assert(decVar > -1);

        lp->addColToRow(rowSmallerThan, decVar, m);

        for (size_t p = 0; p < segment->pl().getCardinality(); ++p) {
          std::stringstream ss;
          ss << "x_(" << segment->pl().getStrRepr()
             << ",l=" << linepair.first.line << ",p<=" << p << ")";

          int first = lp->getVarByName(ss.str());
          assert(first > -1);

          std::stringstream ss2;
          ss2 << "x_(" << segment->pl().getStrRepr()
              << ",l=" << linepair.second.line << ",p<=" << p << ")";

          int second = lp->getVarByName(ss2.str());
          assert(second > -1);

          lp->addColToRow(rowSmallerThan, first, 1);
          lp->addColToRow(rowSmallerThan, second, -1);
        }
      }
    }
  }

  // sum constraint for separation
  for (OptNode* node : g) {
    for (OptEdge* segment : node->getAdjList()) {
      if (segment->getFrom() != node) continue;
      for (LinePair linepair : getLinePairs(segment, true)) {
        std::stringstream rowName;

        int rowDistance1 = 0;
        int rowDistance2 = 0;
        if (_cfg->splittingOpt && segment->pl().getCardinality() > 2) {
          rowName.str("");
          rowName << "sum_distancor1(e=" << segment->pl().getStrRepr()
                  << ",A=" << linepair.first.line
                  << ",B=" << linepair.second.line << ")";

          rowDistance1 = lp->addRow(rowName.str(), 1, shared::optim::UP);

          rowName.str("");
          rowName << "sum_distancor2(e=" << segment->pl().getStrRepr()
                  << ",A=" << linepair.first.line
                  << ",B=" << linepair.second.line << ")";

          rowDistance2 = lp->addRow(rowName.str(), 1, shared::optim::UP);

          std::stringstream sss;
          sss << "x_(" << segment->pl().getStrRepr() << ","
              << linepair.first.line << "<T>" << linepair.second.line << ")";

          int decVarDistance = lp->getVarByName(sss.str());
          assert(decVarDistance > -1);

          lp->addColToRow(rowDistance1, decVarDistance, -static_cast<int>(m));
          lp->addColToRow(rowDistance2, decVarDistance, -static_cast<int>(m));

          for (size_t p = 0; p < segment->pl().getCardinality(); ++p) {
            std::stringstream ss;
            ss << "x_(" << segment->pl().getStrRepr()
               << ",l=" << linepair.first.line << ",p<=" << p << ")";

            int first = lp->getVarByName(ss.str());
            assert(first > -1);

            std::stringstream ss2;
            ss2 << "x_(" << segment->pl().getStrRepr()
                << ",l=" << linepair.second.line << ",p<=" << p << ")";

            int second = lp->getVarByName(ss2.str());
            assert(second > -1);

            lp->addColToRow(rowDistance1, first, 1);
            lp->addColToRow(rowDistance1, second, -1);

            lp->addColToRow(rowDistance2, first, -1);
            lp->addColToRow(rowDistance2, second, 1);
          }
        }
      }
    }
  }

  // crossing constraints
  for (OptNode* node : g) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->getAdjList()) {
      processed.insert(segmentA);

      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA, true)) {
        // iterate over all edges this
        // pair traverses to _TOGETHER_
        // (its possible that there are multiple edges if a line continues
        //  in more then 1 segment)
        for (OptEdge* segmentB : getEdgePartners(node, segmentA, linepair)) {
          if (processed.find(segmentB) != processed.end()) continue;

          // introduce dec var
          std::stringstream ss;
          ss << "x_dec(" << segmentA->pl().getStrRepr() << ","
             << segmentA->pl().getStrRepr() << segmentB->pl().getStrRepr()
             << "," << linepair.first.line << "(" << linepair.first.line->id()
             << ")," << linepair.second.line << "("
             << linepair.second.line->id() << ")," << node << ")";

          int decisionVar = lp->addCol(
              ss.str(), shared::optim::BIN,
              getCrossingPenaltySameSeg(node)
                  // multiply the penalty with the number of collapsed lines!
                  * (linepair.first.relatives.size()) *
                  (linepair.second.relatives.size()));

          int aSmallerBinL1 = 0;
          int aSmallerBinL2 = 0;
          int bSmallerAinL2 = 0;

          std::stringstream aSmBStr;
          aSmBStr << "x_(" << segmentA->pl().getStrRepr() << ","
                  << linepair.first.line << "<" << linepair.second.line << ")";
          aSmallerBinL1 = lp->getVarByName(aSmBStr.str());

          assert(aSmallerBinL1 > -1);

          std::stringstream aBgBStr;
          aBgBStr << "x_(" << segmentB->pl().getStrRepr() << ","
                  << linepair.first.line << "<" << linepair.second.line << ")";
          aSmallerBinL2 = lp->getVarByName(aBgBStr.str());

          assert(aSmallerBinL2 > -1);

          std::stringstream bBgAStr;
          bBgAStr << "x_(" << segmentB->pl().getStrRepr() << ","
                  << linepair.second.line << "<" << linepair.first.line << ")";
          bSmallerAinL2 = lp->getVarByName(bBgAStr.str());

          assert(bSmallerAinL2 > -1);

          std::stringstream rowName;
          rowName << "sum_dec(e1=" << segmentA->pl().getStrRepr()
                  << ",e2=" << segmentB->pl().getStrRepr()
                  << ",A=" << linepair.first.line
                  << ",B=" << linepair.second.line << ",n=" << node << ")";

          int row = lp->addRow(rowName.str(), 0, shared::optim::LO);

          std::stringstream rowName2;
          rowName2 << "sum_dec2(e1=" << segmentA->pl().getStrRepr()
                   << ",e2=" << segmentB->pl().getStrRepr()
                   << ",A=" << linepair.first.line
                   << ",B=" << linepair.second.line << ",n=" << node << ")";

          int row2 = lp->addRow(rowName2.str(), 0, shared::optim::LO);

          bool otherWayA =
              (segmentA->getFrom() != node) ^ segmentA->pl().etgs.front().dir;
          bool otherWayB =
              (segmentB->getFrom() != node) ^ segmentB->pl().etgs.front().dir;

          if (otherWayA ^ otherWayB) {
          } else {
            aSmallerBinL2 = bSmallerAinL2;
          }

          lp->addColToRow(row, aSmallerBinL1, -1);
          lp->addColToRow(row, aSmallerBinL2, 1);
          lp->addColToRow(row, decisionVar, 1);

          lp->addColToRow(row2, aSmallerBinL1, 1);
          lp->addColToRow(row2, aSmallerBinL2, -1);
          lp->addColToRow(row2, decisionVar, 1);
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
            if (segmentA->pl().getCardinality() > 2 &&
                segmentB->pl().getCardinality() > 2) {
              // the interesting case where the line continue together from
              // segment A to segment B and the cardinality of both A and B
              // is > 2 (that is, it is possible in A or B that the two lines
              // won't be together)
              std::stringstream sss;
              sss << "x_decT(" << segmentA->pl().getStrRepr() << ","
                  << segmentA->pl().getStrRepr() << segmentB->pl().getStrRepr()
                  << "," << linepair.first.line << "("
                  << linepair.first.line->id() << ")," << linepair.second.line
                  << "(" << linepair.second.line->id() << ")," << node << ")";

              int decisionVarDist1Change = lp->addCol(
                  sss.str(), shared::optim::BIN, getSplittingPenalty(node));

              int aNearBinL1 = 0;
              int aNearBinL2 = 0;

              std::stringstream aNearBL1Str;
              aNearBL1Str << "x_(" << segmentA->pl().getStrRepr() << ","
                          << linepair.first.line << "<T>"
                          << linepair.second.line << ")";
              aNearBinL1 = lp->getVarByName(aNearBL1Str.str());
              assert(aNearBinL1 > -1);

              std::stringstream aNearBL2Str;
              aNearBL2Str << "x_(" << segmentB->pl().getStrRepr() << ","
                          << linepair.first.line << "<T>"
                          << linepair.second.line << ")";
              aNearBinL2 = lp->getVarByName(aNearBL2Str.str());
              assert(aNearBinL2 > -1);

              std::stringstream rowTName;
              rowTName << "sum_decT(e1=" << segmentA->pl().getStrRepr()
                       << ",e2=" << segmentB->pl().getStrRepr()
                       << ",A=" << linepair.first.line
                       << ",B=" << linepair.second.line << ",n=" << node << ")";

              int rowT = lp->addRow(rowTName.str(), 0, shared::optim::LO);

              std::stringstream rowTName2;
              rowTName2 << "sum_decT2(e1=" << segmentA->pl().getStrRepr()
                        << ",e2=" << segmentB->pl().getStrRepr()
                        << ",A=" << linepair.first.line
                        << ",B=" << linepair.second.line << ",n=" << node
                        << ")";

              int rowT2 = lp->addRow(rowTName2.str(), 0, shared::optim::LO);

              lp->addColToRow(rowT, aNearBinL1, -1);
              lp->addColToRow(rowT, aNearBinL2, 1);
              lp->addColToRow(rowT, decisionVarDist1Change, 1);

              lp->addColToRow(rowT2, aNearBinL1, 1);
              lp->addColToRow(rowT2, aNearBinL2, -1);
              lp->addColToRow(rowT2, decisionVarDist1Change, 1);
            } else if ((segmentA->pl().getCardinality() == 2) ^
                       (segmentB->pl().getCardinality() == 2)) {
              // the trivial case where one of the two segments only has
              // cardinality = 2, so the lines will always be together

              OptEdge* segment =
                  segmentA->pl().getCardinality() != 2 ? segmentA : segmentB;

              std::stringstream aNearBStr;
              aNearBStr << "x_(" << segment->pl().getStrRepr() << ","
                        << linepair.first.line << "<T>" << linepair.second.line
                        << ")";

              lp->setObjCoef(aNearBStr.str(), getSplittingPenalty(node));
            }
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ILPEdgeOrderOptimizer::writeDiffSegConstraintsImpr(
    const std::set<OptNode*>& g, ILPSolver* lp) const {
  // go into nodes and build crossing constraints for adjacent
  for (OptNode* node : g) {
    std::set<OptEdge*> processed;
    for (OptEdge* segmentA : node->getAdjList()) {
      processed.insert(segmentA);
      // iterate over all possible line pairs in this segment
      for (LinePair linepair : getLinePairs(segmentA, true)) {
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
              int testVar = 0;

              if (poscomb.first > poscomb.second) {
                std::stringstream bBgAStr;
                bBgAStr << "x_(" << segmentA->pl().getStrRepr() << ","
                        << linepair.first.line << "<" << linepair.second.line
                        << ")";
                testVar = lp->getVarByName(bBgAStr.str());
              } else {
                std::stringstream bBgAStr;
                bBgAStr << "x_(" << segmentA->pl().getStrRepr() << ","
                        << linepair.second.line << "<" << linepair.first.line
                        << ")";
                testVar = lp->getVarByName(bBgAStr.str());
              }

              assert(testVar);
              std::stringstream ss;
              ss << "dec_sum(" << segmentA->pl().getStrRepr() << ","
                 << segments.first->pl().getStrRepr()
                 << segments.second->pl().getStrRepr() << ","
                 << linepair.first.line << "," << linepair.second.line
                 << "pa=" << poscomb.first << ",pb=" << poscomb.second
                 << ",n=" << node << ")";

              int row = lp->addRow(ss.str(), 0, shared::optim::FIX);

              lp->addColToRow(row, testVar, 1);
              lp->addColToRow(row, decisionVar, -1);

              // one cross is enough...
              break;
            }
          }
        }
      }
    }
  }
}
