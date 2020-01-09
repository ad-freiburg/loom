// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/linegraph/Route.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/Scorer.h"
#include "util/Misc.h"

using namespace transitmapper;
using namespace optim;
using transitmapper::graph::TransitGraph;
using shared::linegraph::Route;
using shared::linegraph::InnerGeometry;
using transitmapper::graph::IDENTITY_PENALTIES;

// _____________________________________________________________________________
double Scorer::getScore() const { return getScore(_g->getConfig()); }

// _____________________________________________________________________________
double Scorer::getScore(const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getScore(n, c);
  }

  return ret;
}
//
// _____________________________________________________________________________
double Scorer::getCrossScore() const { return getCrossScore(_g->getConfig()); }

// _____________________________________________________________________________
double Scorer::getCrossScore(const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getCrossingScore(n, c, _pens);
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getSeparationScore() const {
  return getSeparationScore(_g->getConfig());
}

// _____________________________________________________________________________
double Scorer::getSeparationScore(const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getSeparationScore(n, c, _pens);
  }

  return ret;
}

// _____________________________________________________________________________
size_t Scorer::getNumCrossings() const {
  return getNumCrossings(_g->getConfig());
}

// _____________________________________________________________________________
size_t Scorer::getNumCrossings(const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getNumCrossings(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
size_t Scorer::getNumSeparations() const {
  return getNumSeparations(_g->getConfig());
}

// _____________________________________________________________________________
size_t Scorer::getNumSeparations(const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNds()) {
    ret += getNumSeparations(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getNumPossSolutions() const {
  double ret = 1;

  for (auto n : _g->getNds()) {
    for (auto e : n->getAdjListOut()) {
      ret *= util::factorial(e->pl().getRoutes().size());
    }
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getScore(const shared::linegraph::LineNode* n,
                        const graph::OrderingConfig& cfg) const {
  return getCrossingScore(n, cfg, _pens) + getSeparationScore(n, cfg, _pens);
}

// _____________________________________________________________________________
size_t Scorer::getNumCrossings(const shared::linegraph::LineNode* n,
                               const OrderingConfig& c) const {
  return getCrossingScore(n, c, IDENTITY_PENALTIES);
}

// _____________________________________________________________________________
double Scorer::getCrossingScore(const shared::linegraph::LineNode* n,
                                const OrderingConfig& c,
                                const Penalties& pens) const {
  std::vector<InnerGeometry> igs = _g->getInnerGeometries(n, c, -1);
  size_t ret = 0;

  for (size_t i = 0; i < igs.size(); ++i) {
    for (size_t j = i + 1; j < igs.size(); ++j) {
      const InnerGeometry& iga = igs[i];
      const InnerGeometry& igb = igs[j];

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

  return ret;
}

// _____________________________________________________________________________
size_t Scorer::getNumSeparations(const shared::linegraph::LineNode* n,
                                 const OrderingConfig& c) const {
  return getSeparationScore(n, c, IDENTITY_PENALTIES);
}

// _____________________________________________________________________________
double Scorer::getSeparationScore(const shared::linegraph::LineNode* n,
                                  const OrderingConfig& c,
                                  const Penalties& pens) const {
  size_t ret = 0;
  for (auto nf : n->pl().getMainDirs()) {
    const auto* e = nf.edge;
    std::vector<std::pair<const Route*, const Route*> > curPairs;
    for (size_t i = 0; i < c.find(e)->second.size() - 1; i++) {
      size_t p = c.find(e)->second[i];
      curPairs.push_back(std::pair<const Route*, const Route*>(
          e->pl().routeOccAtPos(p).route,
          e->pl().routeOccAtPos(c.find(e)->second[i + 1]).route));
    }

    for (auto p : curPairs) {
      for (auto nf : n->pl().getMainDirs()) {
        const auto* f = nf.edge;
        if (e == f) continue;

        if (f->pl().hasRoute(p.first) && f->pl().hasRoute(p.second) &&
            n->pl().connOccurs(p.first, e, f) &&
            n->pl().connOccurs(p.second, e, f)) {
          if (abs(int(f->pl().getRoutePosUnder(p.first, c.find(f)->second)) -
                  int(f->pl().getRoutePosUnder(p.second, c.find(f)->second))) >
              1) {
            ret += getSplittingPenalty(n, pens);
          }
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltySameSeg(const shared::linegraph::LineNode* n,
                                      const Penalties& pens) const {
  double ret = 1;
  if (pens.crossAdjPen)
    ret *= n->getAdjListOut().size() + n->getAdjListIn().size();

  if (n->pl().getStops().size() > 0) {
    if (n->getAdjListOut().size() + n->getAdjListIn().size() == 2)
      return pens.inStatCrossPenDegTwo;
    return pens.inStatCrossPenSameSeg * ret;
  }

  return ret * pens.sameSegCrossPen;
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltyDiffSeg(const shared::linegraph::LineNode* n,
                                      const Penalties& pens) const {
  double ret = 1;
  if (pens.crossAdjPen)
    ret *= n->getAdjListOut().size() + n->getAdjListIn().size();

  if (n->pl().getStops().size() > 0) {
    if (n->getAdjListOut().size() + n->getAdjListIn().size() == 2)
      return pens.inStatCrossPenDegTwo;
    return pens.inStatCrossPenDiffSeg * ret;
  }

  return ret * pens.diffSegCrossPen;
}

// _____________________________________________________________________________
int Scorer::getSplittingPenalty(const shared::linegraph::LineNode* n,
                                const Penalties& pens) const {
  double ret = 1;
  if (pens.splitAdjPen)
    ret *= n->getAdjListOut().size() + n->getAdjListIn().size();

  if (n->pl().getStops().size() > 0) {
    if (n->getAdjListOut().size() + n->getAdjListIn().size() == 2)
      return pens.inStatSplitPenDegTwo;
    return pens.inStatSplitPen * ret;
  }

  return ret * pens.splitPen;
}

// _____________________________________________________________________________
double Scorer::getSeparationScore(const shared::linegraph::LineNode* n,
                                  const OrderingConfig& c) const {
  return getSeparationScore(n, c, _pens);
}

// _____________________________________________________________________________
double Scorer::getCrossingScore(const shared::linegraph::LineNode* n,
                                const OrderingConfig& c) const {
  return getCrossingScore(n, c, _pens);
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltySameSeg(
    const shared::linegraph::LineNode* n) const {
  return getCrossingPenaltySameSeg(n, _pens);
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltyDiffSeg(
    const shared::linegraph::LineNode* n) const {
  return getCrossingPenaltyDiffSeg(n, _pens);
}

// _____________________________________________________________________________
int Scorer::getSplittingPenalty(const shared::linegraph::LineNode* n) const {
  return getSplittingPenalty(n, _pens);
}
