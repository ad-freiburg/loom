// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/graph/Penalties.h"
#include "transitmap/optim/OptGraph.h"
#include "transitmap/optim/OptGraphScorer.h"
#include "transitmap/optim/Optimizer.h"

using namespace transitmapper;
using namespace optim;
using transitmapper::graph::TransitGraph;
using transitmapper::graph::Node;
using transitmapper::graph::Edge;
using transitmapper::graph::Route;
using transitmapper::graph::InnerGeometry;
using transitmapper::graph::IDENTITY_PENALTIES;

// _____________________________________________________________________________
double OptGraphScorer::getScore(const std::set<OptNode*>& g,
                                const OptOrderingConfig& c) const {
  double ret = 0;

  for (auto n : g) {
    ret += getScore(n, c);
  }

  return ret;
}

// _____________________________________________________________________________
double OptGraphScorer::getScore(OptNode* n,
                                const OptOrderingConfig& c) const {
  if (!n->pl().node) return 0;
  auto numCrossings = getNumCrossings(n, c);

  // std::cout << "Number of crossings in node " << n;
  // if (n->pl().node && n->pl().node->getStops().size()) std::cout << " (" << n->pl().node->getStops().front().name << ")" << std::endl;

  // std::cout << ": " << numCrossings.first
            // << ", " << numCrossings.second << std::endl;

  return numCrossings.first * _scorer->getCrossingPenaltySameSeg(n->pl().node) + numCrossings.second * _scorer->getCrossingPenaltyDiffSeg(n->pl().node);
}

// _____________________________________________________________________________
std::pair<size_t, size_t> OptGraphScorer::getNumCrossings(
    OptNode* n, const OptOrderingConfig& c) const {
  size_t sameSegCrossings = 0;
  size_t diffSegCrossings = 0;

  std::map<const Route*, std::set<OptEdge*>> proced;

  for (auto ea : n->getAdjList()) {
    auto linePairs = Optimizer::getLinePairs(ea, true);

    for (auto lp : linePairs) {
      // check if pairs continue in same segments
      proced[lp.first].insert(ea);
      proced[lp.second].insert(ea);

      for (auto eb : Optimizer::getEdgePartners(n, ea, lp)) {
        if (proced[lp.first].count(eb) && proced[lp.second].count(eb)) continue;

        proced[lp.first].insert(eb);
        proced[lp.second].insert(eb);

        PosCom posA(
            std::distance(c.at(ea).begin(),
                          std::find(c.at(ea).begin(), c.at(ea).end(), lp.first)),
            std::distance(c.at(eb).begin(),
                          std::find(c.at(eb).begin(), c.at(eb).end(), lp.first)));

        PosCom posB(
            std::distance(c.at(ea).begin(),
                          std::find(c.at(ea).begin(), c.at(ea).end(), lp.second)),
            std::distance(c.at(eb).begin(),
                          std::find(c.at(eb).begin(), c.at(eb).end(), lp.second)));



        PosComPair poses(posA, posB);

        if (Optimizer::crosses(n, ea, eb, poses)) {
          sameSegCrossings++;
        }
      }

      for (auto ebc : Optimizer::getEdgePartnerPairs(n, ea, lp)) {
        PosCom posA(
            std::distance(c.at(ea).begin(),
                          std::find(c.at(ea).begin(), c.at(ea).end(), lp.first)),
            std::distance(
                c.at(ea).begin(),
                std::find(c.at(ea).begin(), c.at(ea).end(), lp.second)));

        if (Optimizer::crosses(n, ea, ebc, posA)) diffSegCrossings++;
      }
    }
  }

  return std::pair<size_t, size_t>(sameSegCrossings, diffSegCrossings);
}
