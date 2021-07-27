// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "loom/optim/Optimizer.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/Penalties.h"

using namespace loom;
using namespace optim;
using shared::linegraph::InnerGeom;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;

// _____________________________________________________________________________
double OptGraphScorer::getSplittingScore(OptGraph* og,
                                         const std::set<OptNode*>& g,
                                         const OptOrderCfg& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getSplittingScore(og, n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingScore(OptGraph* og,
                                        const std::set<OptNode*>& g,
                                        const OptOrderCfg& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getCrossingScore(og, n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingScore(OptGraph* og, OptNode* n,
                                        const OptOrderCfg& c) const {
  if (!n->pl().node) return 0;
  auto numCrossings = getNumCrossings(og, n, c);

  return numCrossings.first * _scorer->getCrossingPenaltySameSeg(n->pl().node) +
         numCrossings.second * _scorer->getCrossingPenaltyDiffSeg(n->pl().node);
}

// _____________________________________________________________________________
double OptGraphScorer::getSplittingScore(OptGraph* og, OptNode* n,
                                         const OptOrderCfg& c) const {
  if (!n->pl().node) return 0;
  return getNumSeparations(og, n, c) *
         _scorer->getSplittingPenalty(n->pl().node);
}

// _____________________________________________________________________________
size_t OptGraphScorer::getNumSeparations(OptGraph* og, OptNode* n,
                                         const OptOrderCfg& c) const {
  size_t seps = 0;

  for (auto ea : n->getAdjList()) {
    auto linePairs = Optimizer::getLinePairs(ea, true);

    for (auto lp : linePairs) {
      for (auto eb : Optimizer::getEdgePartners(n, ea, lp)) {
        int ainA = std::distance(
            c.at(ea).begin(),
            std::find(c.at(ea).begin(), c.at(ea).end(), lp.first.line));
        int ainB = std::distance(
            c.at(eb).begin(),
            std::find(c.at(eb).begin(), c.at(eb).end(), lp.first.line));

        int binA = std::distance(
            c.at(ea).begin(),
            std::find(c.at(ea).begin(), c.at(ea).end(), lp.second.line));
        int binB = std::distance(
            c.at(eb).begin(),
            std::find(c.at(eb).begin(), c.at(eb).end(), lp.second.line));

        if (abs(ainA - binA) == 1 && abs(ainB - binB) != 1) {
          seps++;
        }
      }
    }
  }

  return seps;
}

// _____________________________________________________________________________
std::pair<size_t, size_t> OptGraphScorer::getNumCrossings(
    OptGraph* og, OptNode* n, const OptOrderCfg& c) const {
  size_t sameSegCrossings = 0;
  size_t diffSegCrossings = 0;

  std::map<LinePair, std::set<OptEdge*>> proced;

  for (auto ea : n->getAdjList()) {
    // line pairs are unique because of the second parameter
    // they are always sorted by their pointer value
    auto linePairs = Optimizer::getLinePairs(ea, true);

    for (auto lp : linePairs) {
      // check if pairs continue in same segments

      // mark this line pair as processed on ea - we have checked it
      // into each adjacent edge
      proced[lp].insert(ea);

      for (auto eb : Optimizer::getEdgePartners(n, ea, lp)) {
        // if we have already fully checked the line pairs on this edge,
        // don't count the crossing again - skip.
        if (proced[lp].count(eb)) continue;

        PosCom posA(std::distance(c.at(ea).begin(),
                                  std::find(c.at(ea).begin(), c.at(ea).end(),
                                            lp.first.line)),
                    std::distance(c.at(eb).begin(),
                                  std::find(c.at(eb).begin(), c.at(eb).end(),
                                            lp.first.line)));

        PosCom posB(std::distance(c.at(ea).begin(),
                                  std::find(c.at(ea).begin(), c.at(ea).end(),
                                            lp.second.line)),
                    std::distance(c.at(eb).begin(),
                                  std::find(c.at(eb).begin(), c.at(eb).end(),
                                            lp.second.line)));

        PosComPair poses(posA, posB);

        if (Optimizer::crosses(n, ea, eb, poses)) sameSegCrossings++;
      }

      for (auto ebc : Optimizer::getEdgePartnerPairs(n, ea, lp)) {
        PosCom posA(std::distance(c.at(ea).begin(),
                                  std::find(c.at(ea).begin(), c.at(ea).end(),
                                            lp.first.line)),
                    std::distance(c.at(ea).begin(),
                                  std::find(c.at(ea).begin(), c.at(ea).end(),
                                            lp.second.line)));

        if (Optimizer::crosses(n, ea, ebc, posA)) diffSegCrossings++;
      }
    }
  }

  return std::pair<size_t, size_t>(sameSegCrossings, diffSegCrossings);
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingScore(OptGraph* og, OptEdge* e,
                                        const OptOrderCfg& c) const {
  return getCrossingScore(og, e->getFrom(), c) +
         getCrossingScore(og, e->getTo(), c);
}

// _____________________________________________________________________________
double OptGraphScorer::getSplittingScore(OptGraph* og, OptEdge* e,
                                         const OptOrderCfg& c) const {
  return getSplittingScore(og, e->getFrom(), c) +
         getSplittingScore(og, e->getTo(), c);
}
