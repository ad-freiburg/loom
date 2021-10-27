// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "loom/optim/GreedyOptimizer.h"
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
                                     HierarOrderCfg* hc, size_t depth,
                                     OptResStats& stats) const {
  UNUSED(depth);
  OptOrderCfg cur;

  // fixed order list of optim graph edges
  std::vector<OptEdge*> edges;

  for (auto n : g)
    for (auto e : n->getAdjList())
      if (n == e->getFrom() && e->pl().getCardinality() > 1) edges.push_back(e);

  if (_randomStart) {
    // this is the starting ordering, which is random
    initialConfig(g, &cur, false);
  } else {
    // take the greedy optimized ordering as a starting point
    GreedyOptimizer greedy(_cfg, _scorer.getPens(), true);
    greedy.getFlatConfig(g, &cur);
  }

  size_t iters = 0;
  size_t last = 0;

  while (true) {
    iters++;

    bool found = false;
    double bestChange = 0;
    OptEdge* bestEdge;
    std::vector<const Line*> bestOrder;

    for (size_t i = 0; i < edges.size(); i++) {
      double oldScore = getScore(og, edges[i], cur);

      for (size_t p1 = 0; p1 < cur[edges[i]].size(); p1++) {
        for (size_t p2 = p1; p2 < cur[edges[i]].size(); p2++) {
          // switch p1 and p2
          auto tmp = cur[edges[i]][p1];
          cur[edges[i]][p1] = cur[edges[i]][p2];
          cur[edges[i]][p2] = tmp;

          double s = getScore(og, edges[i], cur);
          if (s < oldScore && oldScore - s > bestChange) {
            found = true;
            bestChange = oldScore - s;
            bestEdge = edges[i];
            bestOrder = cur[edges[i]];
          }

          // switch back
          tmp = cur[edges[i]][p1];
          cur[edges[i]][p1] = cur[edges[i]][p2];
          cur[edges[i]][p2] = tmp;
        }
      }
    }

    if (!found) break;

    cur[bestEdge] = bestOrder;
  }

  writeHierarch(&cur, hc);
  return iters;
}

// _____________________________________________________________________________
double HillClimbOptimizer::getScore(OptGraph* og, OptEdge* e,
                                    OptOrderCfg& cur) const {
  if (_optScorer.optimizeSep()) return _optScorer.getTotalScore(e, cur);
  return _optScorer.getCrossingScore(e, cur);
}
