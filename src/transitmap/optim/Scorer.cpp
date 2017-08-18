// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/Scorer.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/graph/Node.h"
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/Route.h"
#include "transitmap/graph/Penalties.h"
#include "util/Misc.h"

using namespace transitmapper;
using namespace optim;
using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using transitmapper::graph::Route;
using transitmapper::graph::InnerGeometry;
using transitmapper::graph::IDENTITY_PENALTIES;

// _____________________________________________________________________________
double Scorer::getScore(const Penalties& pens) const {
  return getScore(pens, _g->getConfig());
}

// _____________________________________________________________________________
double Scorer::getScore(const Penalties& pens, const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNodes()) {
    ret += getScore(n, pens, true, true, c);
  }

  return ret;
}
//
// _____________________________________________________________________________
double Scorer::getCrossScore(const Penalties& pens) const {
  return getCrossScore(pens, _g->getConfig());
}

// _____________________________________________________________________________
double Scorer::getCrossScore(const Penalties& pens, const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNodes()) {
    ret += getCrossingScore(n, c, pens, true);
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getSeparationScore(const Penalties& pens) const {
  return getSeparationScore(pens, _g->getConfig());
}

// _____________________________________________________________________________
double Scorer::getSeparationScore(const Penalties& pens,
                                        const OrderingConfig& c) const {
  double ret = 0;

  for (auto n : _g->getNodes()) {
    ret += getSeparationScore(n, c, pens, true);
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

  for (auto n : _g->getNodes()) {
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

  for (auto n : _g->getNodes()) {
    ret += getNumSeparations(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getNumPossSolutions() const {
  double ret = 1;

  for (auto n : _g->getNodes()) {
    for (auto e : n->getAdjListOut()) {
      ret *= util::factorial(e->getCardinality());
    }
  }

  return ret;
}

// _____________________________________________________________________________
double Scorer::getScore(const Node* n, const Penalties& pens,
                      bool crossAdjPen, bool splitAdjPen,
                      const graph::OrderingConfig& cfg) const {
  return getCrossingScore(n, cfg, pens, crossAdjPen) +
         getSeparationScore(n, cfg, pens, splitAdjPen);
}

// _____________________________________________________________________________
size_t Scorer::getNumCrossings(const Node* n, const OrderingConfig& c) const {
  return getCrossingScore(n, c, IDENTITY_PENALTIES, false);
}

// _____________________________________________________________________________
double Scorer::getCrossingScore(const Node* n, const OrderingConfig& c, const Penalties& pens,
                              bool adjpen) const {
  std::vector<InnerGeometry> igs = n->getInnerGeometries(c, -1);
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

      if (util::geo::intersects(
              iga.geom.getLine().front(), iga.geom.getLine().back(),
              igb.geom.getLine().front(), igb.geom.getLine().back()) ||
          util::geo::dist(iga.geom.getLine(), igb.geom.getLine()) < 1) {
        if (sameSeg) {
          ret += 1 * getCrossingPenaltySameSeg(n, pens, adjpen);
        } else {
          ret += 1 * getCrossingPenaltyDiffSeg(n, pens, adjpen);
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
size_t Scorer::getNumSeparations(const Node* n, const OrderingConfig& c) const {
  return getSeparationScore(n, c, IDENTITY_PENALTIES, false);
}

// _____________________________________________________________________________
size_t Scorer::getSeparationScore(const Node* n, const OrderingConfig& c, const Penalties& pens, bool adjpen) const {
  size_t ret = 0;
  for (auto nf : n->getMainDirs()) {
    const Edge* e = nf.edge;
    std::vector<std::pair<const Route*, const Route*> > curPairs;
    for (size_t i = 0; i < c.find(e)->second.size() - 1; i++) {
      size_t p = c.find(e)->second[i];
      curPairs.push_back(std::pair<const Route*, const Route*>(
          e->getTripsUnordered().at(p).route,
          e->getTripsUnordered().at(c.find(e)->second[i + 1]).route));
    }

    for (auto p : curPairs) {
      for (auto nf : n->getMainDirs()) {
        const Edge* f = nf.edge;
        if (e == f) continue;

        if (f->containsRoute(p.first) && f->containsRoute(p.second) &&
            n->connOccurs(p.first, e, f) && n->connOccurs(p.second, e, f)) {
          if (abs(int(f->getRouteOccWithPosUnder(p.first, c.find(f)->second)
                          .second) -
                  int(f->getRouteOccWithPosUnder(p.second, c.find(f)->second)
                          .second)) > 1) {
            ret += 1 * getSplittingPenalty(n, pens, adjpen);
          }
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltySameSeg(const Node* n, const Penalties& pens, bool adjpen) const {
  double ret = 1;
  if (adjpen) ret *= n->getAdjListIn().size() + n->getAdjListIn().size();

  if (n->getStops().size() > 0) {
    if (n->getAdjListIn().size() + n->getAdjListIn().size() == 2)
      return pens.inStatCrossPenDegTwo;
    return pens.inStatCrossPenSameSeg * ret;
  }

  return ret * pens.sameSegCrossPen;
}

// _____________________________________________________________________________
int Scorer::getCrossingPenaltyDiffSeg(const Node* n, const Penalties& pens, bool adjpen) const {
  double ret = 1;
  if (adjpen) ret *= n->getAdjListIn().size() + n->getAdjListIn().size();

  if (n->getStops().size() > 0) {
    if (n->getAdjListIn().size() + n->getAdjListIn().size() == 2)
      return pens.inStatCrossPenDegTwo;
    return pens.inStatCrossPenDiffSeg * ret;
  }

  return ret * pens.diffSegCrossPen;
}


// _____________________________________________________________________________
int Scorer::getSplittingPenalty(const Node* n, const Penalties& pens, bool adjpen) const {
  double ret = 1;
  if (adjpen) ret *= n->getAdjListIn().size() + n->getAdjListIn().size();

  if (n->getStops().size() > 0) {
    if (n->getAdjListIn().size() + n->getAdjListIn().size() == 2)
      return pens.inStatSplitPenDegTwo;
    return pens.inStatSplitPen * ret;
  }

  return ret * pens.splitPen;
}