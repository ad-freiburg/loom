// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <map>
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
  std::pair<std::pair<size_t, size_t>, size_t> ret = {{0, 0}, 0};
  for (auto ea : n->getAdjList()) {
    auto cur = getNumCrossSeps(n, ea, c);
    ret.first.first += cur.first.first;
    ret.second += cur.second;
  }

  if (n->getDeg() > 2) {
    // diff seg crossings
    for (auto ea : n->getAdjList()) {
      size_t cur = getNumCrossDiffSeg(n, ea, c);
      ret.first.second += cur;
    }

    ret.first.second -= ret.first.first;
  }

  // same seg crossings are counted twice!
  ret.first.first /= 2;

  return ret;
}

// _____________________________________________________________________________
size_t OptGraphScorer::getNumCrossDiffSeg(OptNode* n, OptEdge* ea,
                                          const OptOrderCfg& c) const {
  std::map<const Line*, size_t> ordering;

  bool revA = (ea->getFrom() != n) ^ ea->pl().lnEdgParts.front().dir;

  const auto& cea = c.at(ea);

  for (size_t i = 0; i < cea.size(); i++) {
    const auto& l = cea[i];
    ordering[l] = revA ? cea.size() - 1 - i : i;
  }

  std::vector<size_t> relOrderCross;

  for (const auto& eb : OptGraph::clockwEdges(ea, n)) {
    const auto& ceb = c.at(eb);
    bool revB = (eb->getFrom() != n) ^ eb->pl().lnEdgParts.front().dir;

    for (size_t i = 0; i < ceb.size(); i++) {
      const auto& thisLine = ceb[!revB ? ceb.size() - 1 - i : i];

      auto otherIt = ordering.find(thisLine);
      if (otherIt == ordering.end()) continue;

      const auto* eaLo = ea->pl().getLineOcc(otherIt->first);
      const auto* ebLo = eb->pl().getLineOcc(thisLine);

      assert(eaLo);
      assert(ebLo);

      if ((eaLo->dir == 0 || ebLo->dir == 0 ||
           (eaLo->dir == n->pl().node && ebLo->dir != n->pl().node) ||
           (eaLo->dir != n->pl().node && ebLo->dir == n->pl().node)) &&
          (n->pl().node->pl().connOccurs(eaLo->line, OptGraph::getAdjEdg(ea, n),
                                         OptGraph::getAdjEdg(eb, n)))) {
        // connection occurs, consider for crossings
        relOrderCross.push_back(otherIt->second);
      }
    }
  }

  return util::inversions(relOrderCross);
}

// _____________________________________________________________________________
std::pair<std::pair<size_t, size_t>, size_t> OptGraphScorer::getNumCrossSeps(
    OptNode* n, OptEdge* ea, OptEdge* eb, const OptOrderCfg& c) const {
  std::pair<std::pair<size_t, size_t>, size_t> ret{{0, 0}, 0};

  std::map<const Line*, size_t> ordering;

  bool revA = (ea->getFrom() != n) ^ ea->pl().lnEdgParts.front().dir;
  bool revB = (eb->getFrom() != n) ^ eb->pl().lnEdgParts.front().dir;

  bool rev = !(revA ^ revB);

  const auto& cea = c.at(ea);
  const auto& ceb = c.at(eb);

  for (size_t i = 0; i < cea.size(); i++) {
    ordering[cea[i]] = rev ? cea.size() - 1 - i : i;
  }

  std::vector<size_t> relOrderCross, relOrderSep;

  for (size_t i = 0; i < ceb.size(); i++) {
    const auto& thisLine = ceb[i];

    auto otherIt = ordering.find(thisLine);
    if (otherIt == ordering.end()) {
      // insert a placeholder for separations, otherwise ignore
      relOrderSep.push_back(std::numeric_limits<size_t>::max());
      continue;
    }

    const auto* eaLo = ea->pl().getLineOcc(otherIt->first);
    const auto* ebLo = eb->pl().getLineOcc(thisLine);

    assert(eaLo);
    assert(ebLo);

    if ((eaLo->dir == 0 || ebLo->dir == 0 ||
         (eaLo->dir == n->pl().node && ebLo->dir != n->pl().node) ||
         (eaLo->dir != n->pl().node && ebLo->dir == n->pl().node)) &&
        (n->pl().node->pl().connOccurs(eaLo->line, OptGraph::getAdjEdg(ea, n),
                                       OptGraph::getAdjEdg(eb, n)))) {
      // connection occurs, consider for crossings
      relOrderCross.push_back(otherIt->second);
      relOrderSep.push_back(otherIt->second);
    } else {
      // otherwise insert a placeholder
      relOrderSep.push_back(std::numeric_limits<size_t>::max());
    }
  }

  size_t seps = 0;

  // count separations
  for (size_t i = 1; i < relOrderSep.size(); i++) {
    if (relOrderSep[i - 1] < std::numeric_limits<size_t>::max() &&
        relOrderSep[i] < std::numeric_limits<size_t>::max()) {
      if (relOrderSep[i] > relOrderSep[i - 1] &&
          relOrderSep[i] - relOrderSep[i - 1] > 1)
        seps++;
      if (relOrderSep[i] < relOrderSep[i - 1] &&
          relOrderSep[i - 1] - relOrderSep[i] > 1)
        seps++;
    }
  }

  ret.first.first = util::inversions(relOrderCross);
  ret.second = seps;

  return ret;
}

// _____________________________________________________________________________
std::pair<std::pair<size_t, size_t>, size_t> OptGraphScorer::getNumCrossSeps(
    OptNode* n, OptEdge* ea, const OptOrderCfg& c) const {
  std::pair<std::pair<size_t, size_t>, size_t> ret = {{0, 0}, 0};
  for (auto eb : n->getAdjList()) {
    if (eb == ea) continue;
    auto cur = getNumCrossSeps(n, ea, eb, c);
    ret.first.first += cur.first.first;
    ret.first.second += cur.first.second;
    ret.second += cur.second;
  }

  return ret;
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
