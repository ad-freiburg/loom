// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "loom/optim/OptGraph.h"
#include "loom/optim/OptGraphScorer.h"
#include "loom/optim/Optimizer.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/Penalties.h"

using loom::optim::OptGraphScorer;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;

// _____________________________________________________________________________
std::pair<size_t, size_t> OptGraphScorer::getNumCrossings(
    const OptGraph* g, const OptOrderCfg& c) const {
  return getNumCrossings(g->getNds(), c);
}

// _____________________________________________________________________________
std::pair<size_t, size_t> OptGraphScorer::getNumCrossings(
    const std::set<OptNode*>& g, const OptOrderCfg& c) const {
  size_t sameSegCrossings = 0;
  size_t diffSegCrossings = 0;

  for (auto n : g) {
    auto crossings = getNumCrossings(n, c);
    sameSegCrossings += crossings.first;
    diffSegCrossings += crossings.second;
  }

  return {sameSegCrossings, diffSegCrossings};
}

// _____________________________________________________________________________
size_t OptGraphScorer::getNumSeparations(const OptGraph* g,
                                         const OptOrderCfg& c) const {
  return getNumSeparations(g->getNds(), c);
}

// _____________________________________________________________________________
size_t OptGraphScorer::getNumSeparations(const std::set<OptNode*>& g,
                                         const OptOrderCfg& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getNumSeparations(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double OptGraphScorer::getSeparationScore(const OptGraph* g,
                                          const OptOrderCfg& c) const {
  return getSeparationScore(g->getNds(), c);
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingScore(const OptGraph* g,
                                        const OptOrderCfg& c) const {
  return getCrossingScore(g->getNds(), c);
}

// _____________________________________________________________________________
double OptGraphScorer::getTotalScore(const OptGraph* g,
                                     const OptOrderCfg& c) const {
  return getTotalScore(g->getNds(), c);
}

// _____________________________________________________________________________
double OptGraphScorer::getTotalScore(const std::set<OptNode*>& g,
                                     const OptOrderCfg& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getTotalScore(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double OptGraphScorer::getTotalScore(OptEdge* e, const OptOrderCfg& c) const {
  return getTotalScore(e->getFrom(), c) + getTotalScore(e->getTo(), c);
}

// _____________________________________________________________________________
double OptGraphScorer::getTotalScore(OptNode* n, const OptOrderCfg& c) const {
  if (!n->pl().node) return 0;
  auto num = getNumCrossSeps(n, c);

  return num.first.first * getCrossingPenSameSeg(n) +
         num.first.second * getCrossingPenDiffSeg(n) +
         num.second * getSeparationPen(n);
}

// _____________________________________________________________________________
double OptGraphScorer::getSeparationScore(const std::set<OptNode*>& g,
                                          const OptOrderCfg& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getSeparationScore(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingScore(const std::set<OptNode*>& g,
                                        const OptOrderCfg& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getCrossingScore(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingScore(OptNode* n,
                                        const OptOrderCfg& c) const {
  if (!n->pl().node) return 0;
  auto numCrossings = getNumCrossings(n, c);

  return numCrossings.first * getCrossingPenSameSeg(n) +
         numCrossings.second * getCrossingPenDiffSeg(n);
}

// _____________________________________________________________________________
double OptGraphScorer::getSeparationScore(OptNode* n,
                                          const OptOrderCfg& c) const {
  if (!n->pl().node) return 0;
  return getNumSeparations(n, c) * getSeparationPen(n);
}

// _____________________________________________________________________________
size_t OptGraphScorer::getNumSeparations(OptNode* n,
                                         const OptOrderCfg& c) const {
  return getNumCrossSeps(n, c).second;
}

// _____________________________________________________________________________
std::pair<size_t, size_t> OptGraphScorer::getNumCrossings(
    OptNode* n, const OptOrderCfg& c) const {
  return getNumCrossSeps(n, c).first;
}

// _____________________________________________________________________________
std::pair<std::pair<size_t, size_t>, size_t> OptGraphScorer::getNumCrossSeps(
    OptNode* n, const OptOrderCfg& c) const {
  if (n->getDeg() == 1) return {{0, 0}, {0}};

  size_t sameSegCrossings = 0;
  size_t diffSegCrossings = 0;
  size_t seps = 0;

  for (auto ea : n->getAdjList()) {
    // line pairs are unique because of the second parameter
    // they are always sorted by their pointer value
    const auto& linePairs = Optimizer::getLinePairs(ea, true);
    const auto& cea = c.at(ea);

    for (const auto& lp : linePairs) {
      // check if pairs continue in same segments

      size_t peaA =
          std::find(cea.begin(), cea.end(), lp.first.line) - cea.begin();
      size_t peaB =
          std::find(cea.begin(), cea.end(), lp.second.line) - cea.begin();

      for (const auto& eb : Optimizer::getEdgePartners(n, ea, lp)) {
        // if we have already fully checked the line pairs on this edge,
        // don't count the crossing again - skip.

        const auto& ceb = c.at(eb);

        PosCom posA(peaA, std::find(ceb.begin(), ceb.end(), lp.first.line) -
                              ceb.begin());

        PosCom posB(peaB, std::find(ceb.begin(), ceb.end(), lp.second.line) -
                              ceb.begin());

        PosComPair poses(posA, posB);

        if (Optimizer::crosses(n, ea, eb, poses)) sameSegCrossings++;
        if (Optimizer::separates(poses)) seps++;
      }

      PosCom posA(peaA, peaB);

      for (const auto& ebc : Optimizer::getEdgePartnerPairs(n, ea, lp)) {
        if (Optimizer::crosses(n, ea, ebc, posA)) diffSegCrossings++;
      }
    }
  }

  return {{sameSegCrossings / 2, diffSegCrossings}, seps / 2};
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingScore(OptEdge* e,
                                        const OptOrderCfg& c) const {
  return getCrossingScore(e->getFrom(), c) + getCrossingScore(e->getTo(), c);
}

// _____________________________________________________________________________
double OptGraphScorer::getSeparationScore(OptEdge* e,
                                          const OptOrderCfg& c) const {
  return getSeparationScore(e->getFrom(), c) +
         getSeparationScore(e->getTo(), c);
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingPenSameSeg(const OptNode* n) const {
  if (n->getDeg() == 1) return 0;
  double ret = 1;
  if (_pens.crossAdjPen) ret *= n->pl().node->getDeg();

  if (n->pl().node->pl().stops().size() > 0) {
    if (n->pl().node->getDeg() == 2) return _pens.inStatCrossPenDegTwo;
    return _pens.inStatCrossPenSameSeg * ret;
  }

  return ret * _pens.sameSegCrossPen;
}

// _____________________________________________________________________________
double OptGraphScorer::getCrossingPenDiffSeg(const OptNode* n) const {
  if (n->getDeg() == 1) return 0;
  double ret = 1;
  if (_pens.crossAdjPen) ret *= n->pl().node->getDeg();

  if (n->pl().node->pl().stops().size() > 0) {
    if (n->pl().node->getDeg() == 2) return _pens.inStatCrossPenDegTwo;
    return _pens.inStatCrossPenDiffSeg * ret;
  }

  return ret * _pens.diffSegCrossPen;
}

// _____________________________________________________________________________
double OptGraphScorer::getSeparationPen(const OptNode* n) const {
  if (n->getDeg() == 1) return 0;
  double ret = 1;
  if (_pens.splitAdjPen) ret *= n->pl().node->getDeg();

  if (n->pl().node->pl().stops().size() > 0) {
    if (n->pl().node->getDeg() == 2) return _pens.inStatSplitPenDegTwo;
    return _pens.inStatSplitPen * ret;
  }

  return ret * _pens.splitPen;
}

// _____________________________________________________________________________
bool OptGraphScorer::optimizeSep() const {
  return _pens.inStatSplitPenDegTwo > 0 || _pens.inStatSplitPen > 0 ||
         _pens.splitPen > 0;
}
