// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "transitmap/optim/HillClimbOptimizer.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::HillClimbOptimizer;

// _____________________________________________________________________________
int HillClimbOptimizer::optimizeComp(const std::set<OptNode*>& g,
                                  HierarchOrderingConfig* hc) const {
  OptOrderingConfig cur, null;

  // fixed order list of optim graph edges
  std::vector<OptEdge*> edges;

  for (auto n : g)
    for (auto e : n->getAdjList())
      if (n == e->getFrom() && e->pl().getCardinality() > 1) edges.push_back(e);

  // this guarantees that all the orderings are sorted!
  initialConfig(g, &null, true);

  // this is the starting ordering, which is random
  initialConfig(g, &cur, false);

  size_t iters = 0;
  size_t last = 0;

  while (true) {
    iters++;

    if (iters - last == 10000) {
      LOG(DEBUG) << "@ " << iters;
      last = iters;
    }

    bool found = false;
    double bestChange = 0;
    OptEdge* bestEdge;
    std::vector<const Route*> bestOrder;

    for (size_t i = 0; i < edges.size(); i++) {
      double oldScore = getScore(edges[i], cur);
      auto old = cur[edges[i]];
      cur[edges[i]] = null[edges[i]];

      do {
        double s = getScore(edges[i], cur);
        if (s < oldScore && oldScore - s > bestChange) {
          found = true;
          bestChange = oldScore - s;
          bestEdge = edges[i];
          bestOrder = cur[edges[i]];
        }
      } while (std::next_permutation(cur[edges[i]].begin(), cur[edges[i]].end()));

      // reverting
      cur[edges[i]] = old;
    }

    double curScore = _optScorer.getCrossingScore(g, cur);
    if (_cfg->splittingOpt) curScore += _optScorer.getSplittingScore(g, cur);

    if (!found) {
      LOG(DEBUG) << "Local optimum found after " << iters << ", target=" << curScore;
      break;
    }

    cur[bestEdge] = bestOrder;

    LOG(DEBUG) << "At round " << iters << ", target=" << curScore;
  }

  writeHierarch(&cur, hc);
  return iters;
}

// _____________________________________________________________________________
double HillClimbOptimizer::getScore(OptEdge* e, OptOrderingConfig& cur) const {
  double curScore = _optScorer.getCrossingScore(e, cur);
  if (_cfg->splittingOpt) curScore += _optScorer.getSplittingScore(e, cur);
  return curScore;
}
