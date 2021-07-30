// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "loom/optim/Scorer.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/Penalties.h"
#include "shared/rendergraph/RenderGraph.h"
#include "util/Misc.h"

using loom::optim::Scorer;
using shared::rendergraph::InnerGeom;
using shared::linegraph::Line;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::rendergraph::IDENTITY_PENALTIES;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::Penalties;
using shared::rendergraph::RenderGraph;

// _____________________________________________________________________________
double Scorer::getScore() const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getScore(n);
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getCrossScore() const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getCrossingScore(n, _pens);
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getSeparationScore() const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getSeparationScore(n, _pens);
  }

  return ret;
}

// _____________________________________________________________________________
size_t Scorer::getNumCrossings() const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getNumCrossings(n);
  }

  return ret;
}

// _____________________________________________________________________________
size_t Scorer::getNumSeparations() const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getNumSeparations(n);
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getNumPossSolutions() const {
  double ret = 1;

  for (auto n : _g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret *= util::factorial(e->pl().getLines().size());
    }
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getScore(const LineNode* n) const {
  return getCrossingScore(n, _pens) + getSeparationScore(n, _pens);
}

// _____________________________________________________________________________
size_t Scorer::getNumCrossings(const LineNode* n) const {
  return getCrossingScore(n, IDENTITY_PENALTIES);
}

// _____________________________________________________________________________
double Scorer::getCrossingScore(const LineNode* n,
                                const Penalties& pens) const {
  std::vector<InnerGeom> igs = _g->innerGeoms(n, -1);
  size_t ret = 0;

  for (size_t i = 0; i < igs.size(); ++i) {
    for (size_t j = i + 1; j < igs.size(); ++j) {
      const InnerGeom& iga = igs[i];
      const InnerGeom& igb = igs[j];

      if (iga.from.front == 0 || iga.to.front == 0 || igb.from.front == 0 ||
          igb.to.front == 0)
        continue;

      if (iga.from.front == igb.from.front && iga.slotFrom == igb.slotFrom)
        continue;
      if (iga.from.front == igb.to.front && iga.slotFrom == igb.slotTo)
        continue;
      if (iga.to.front == igb.to.front && iga.slotTo == igb.slotTo) continue;
      if (iga.to.front == igb.from.front && iga.slotTo == igb.slotFrom)
        continue;

      bool sameSeg =
          (iga.from.front == igb.from.front && iga.to.front == igb.to.front) ||
          (iga.to.front == igb.from.front && iga.from.front == igb.to.front);

      bool unavoidable =
          iga.from.front != igb.from.front && iga.from.front != igb.to.front &&
          iga.to.front != igb.from.front && iga.to.front != igb.to.front;

      bool iSectGeo =
          util::geo::intersects(
              iga.geom.getLine().front(), iga.geom.getLine().back(),
              igb.geom.getLine().front(), igb.geom.getLine().back()) ||
          util::geo::dist(iga.geom.getLine(), igb.geom.getLine()) < 1;

      if (iSectGeo) {
        if (sameSeg) {
          ret += getCrossingPenaltySameSeg(n, pens);
        } else if (!unavoidable) {
          ret += getCrossingPenaltyDiffSeg(n, pens);
        }
      }
    }
  }


  size_t ret2 = 0;

  return ret;
}

// _____________________________________________________________________________
size_t Scorer::getNumSeparations(const LineNode* n) const {
  return getSeparationScore(n, IDENTITY_PENALTIES);
}

// _____________________________________________________________________________
double Scorer::getSeparationScore(const LineNode* n,
                                  const Penalties& pens) const {
  size_t ret = 0;
  for (auto e : n->getAdjList()) {
    std::vector<std::pair<size_t, size_t> > curPairs;

    for (size_t p = 0; p < e->pl().getLines().size() - 1; p++) {
      for (auto f : n->getAdjList()) {
        if (e == f) continue;

        auto lnA = e->pl().lineOccAtPos(p).line;
        auto lnB = e->pl().lineOccAtPos(p + 1).line;

        if (LineGraph::lineCtd(e, f, lnA) && LineGraph::lineCtd(e, f, lnB)) {
          if (abs(int(f->pl().linePos(lnA)) - int(f->pl().linePos(lnB))) > 1) {
            ret += getSplittingPenalty(n, pens);
          }
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltySameSeg(const LineNode* n,
                                      const Penalties& pens) const {
  double ret = 1;
  if (pens.crossAdjPen) ret *= n->getDeg();

  if (n->pl().stops().size() > 0) {
    if (n->getDeg() == 2) return pens.inStatCrossPenDegTwo;
    return pens.inStatCrossPenSameSeg * ret;
  }

  return ret * pens.sameSegCrossPen;
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltyDiffSeg(const LineNode* n,
                                      const Penalties& pens) const {
  double ret = 1;
  if (pens.crossAdjPen) ret *= n->getDeg();

  if (n->pl().stops().size() > 0) {
    if (n->getDeg() == 2) return pens.inStatCrossPenDegTwo;
    return pens.inStatCrossPenDiffSeg * ret;
  }

  return ret * pens.diffSegCrossPen;
}

// _____________________________________________________________________________
int Scorer::getSplittingPenalty(const LineNode* n,
                                const Penalties& pens) const {
  double ret = 1;
  if (pens.splitAdjPen) ret *= n->getDeg();

  if (n->pl().stops().size() > 0) {
    if (n->getDeg() == 2) return pens.inStatSplitPenDegTwo;
    return pens.inStatSplitPen * ret;
  }

  return ret * pens.splitPen;
}

// _____________________________________________________________________________
double Scorer::getSeparationScore(const LineNode* n) const {
  return getSeparationScore(n, _pens);
}

// _____________________________________________________________________________
double Scorer::getCrossingScore(const LineNode* n) const {
  return getCrossingScore(n, _pens);
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltySameSeg(const LineNode* n) const {
  return getCrossingPenaltySameSeg(n, _pens);
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltyDiffSeg(const LineNode* n) const {
  return getCrossingPenaltyDiffSeg(n, _pens);
}

// _____________________________________________________________________________
int Scorer::getSplittingPenalty(const LineNode* n) const {
  return getSplittingPenalty(n, _pens);
}
