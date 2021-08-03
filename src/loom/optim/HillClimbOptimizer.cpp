// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "loom/optim/HillClimbOptimizer.h"
#include "shared/linegraph/Line.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using loom::optim::HillClimbOptimizer;
using shared::linegraph::Line;
using shared::rendergraph::HierarOrderCfg;

// _____________________________________________________________________________
int HillClimbOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                     HierarOrderCfg* hc, size_t depth) const {
  UNUSED(depth);
  OptOrderCfg cur, null;

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
      LOGTO(INFO, std::cerr) << "@ " << iters;
      last = iters;
    }

    bool found = false;
    double bestChange = 0;
    OptEdge* bestEdge;
    std::vector<const Line*> bestOrder;

    for (size_t i = 0; i < edges.size(); i++) {
      double oldScore = getScore(og, edges[i], cur);
      auto old = cur[edges[i]];
      cur[edges[i]] = null[edges[i]];

      int c = util::factorial(edges[i]->pl().getCardinality());

      if (c > 300000) {
        // heuristic: is the neighboordhood is too large, simply take
        // 1000 random choices
        for (size_t k = 0; k < 1000; k++) {
          int steps = rand() % c;
          for (int j = 0; j < steps; j++)
            std::next_permutation(cur[edges[i]].begin(), cur[edges[i]].end());

          double s = getScore(og, edges[i], cur);
          if (s < oldScore && oldScore - s > bestChange) {
            found = true;
            bestChange = oldScore - s;
            bestEdge = edges[i];
            bestOrder = cur[edges[i]];
          }
        }
      } else {
        do {
          double s = getScore(og, edges[i], cur);
          if (s < oldScore && oldScore - s > bestChange) {
            found = true;
            bestChange = oldScore - s;
            bestEdge = edges[i];
            bestOrder = cur[edges[i]];
          }
        } while (
            std::next_permutation(cur[edges[i]].begin(), cur[edges[i]].end()));
      }

      // reverting
      cur[edges[i]] = old;
    }

    double curScore = _optScorer.getCrossingScore(g, cur);
    if (_cfg->splittingOpt) curScore += _optScorer.getSplittingScore(g, cur);

    if (!found) {
      LOGTO(INFO, std::cerr)
          << "Local optimum found after " << iters << ", target=" << curScore;
      break;
    }

    cur[bestEdge] = bestOrder;

    LOGTO(DEBUG, std::cerr) << "At round " << iters << ", target=" << curScore;
  }

  writeHierarch(&cur, hc);
  return iters;
}

// _____________________________________________________________________________
double HillClimbOptimizer::getScore(OptGraph* og, OptEdge* e,
                                    OptOrderCfg& cur) const {
  if (_cfg->splittingOpt) return _optScorer.getTotalScore(e, cur);
  return _optScorer.getCrossingScore(e, cur);
}
