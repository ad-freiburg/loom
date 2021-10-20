// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "loom/optim/GreedyOptimizer.h"
#include "loom/optim/SimulatedAnnealingOptimizer.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using loom::optim::SimulatedAnnealingOptimizer;
using shared::linegraph::Line;
using shared::rendergraph::HierarOrderCfg;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::RenderGraph;

// _____________________________________________________________________________
int SimulatedAnnealingOptimizer::optimizeComp(OptGraph* og,
                                              const std::set<OptNode*>& g,
                                              HierarOrderCfg* hc, size_t depth,
                                              OptResStats& stats) const {
  OptOrderCfg cur;

  // fixed order list of optim graph edges
  std::vector<OptEdge*> edges;

  for (auto n : g)
    for (auto e : n->getAdjList())
      if (n == e->getFrom()) edges.push_back(e);

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

  size_t k = 0;

  size_t ABORT_AFTER_UNCH = 50;

  while (true) {
    iters++;

    double temp = 1000.0 / iters;

    if (iters - last == 10000) {
      LOGTO(DEBUG, std::cerr) << "@ " << iters << ", temp = " << temp
                              << ", last change was at " << k << " iters.";
      last = iters;
    }

    for (size_t i = 0; i < edges.size(); i++) {
      double oldScore = getScore(og, edges[i], cur);

      for (size_t p1 = 0; p1 < cur[edges[i]].size(); p1++) {
        for (size_t p2 = p1; p2 < cur[edges[i]].size(); p2++) {
          // switch p1 and p2
          auto tmp = cur[edges[i]][p1];
          cur[edges[i]][p1] = cur[edges[i]][p2];
          cur[edges[i]][p2] = tmp;

          double s = getScore(og, edges[i], cur);

          double r = rand() / (RAND_MAX + 1.0);
          double e = exp(-(1.0 * (s - oldScore)) / temp);

          if (s < oldScore) {
            // found a better solution, keep it, update score
            oldScore = s;
            k = iters;
          } else if (s != oldScore && e > r) {
            // keep solution, despite not bringing any local gain, update score
            oldScore = s;
            k = iters;
          } else {
            // switch back
            tmp = cur[edges[i]][p1];
            cur[edges[i]][p1] = cur[edges[i]][p2];
            cur[edges[i]][p2] = tmp;
          }
        }
      }
    }

    if (iters - k > ABORT_AFTER_UNCH) break;
  }

  double curScore = 0;
  if (_optScorer.optimizeSep())
    curScore += _optScorer.getTotalScore(g, cur);
  else
    curScore = _optScorer.getCrossingScore(g, cur);

  LOGTO(INFO, std::cerr) << "Stopped after " << iters
                         << " iterations. Final target = " << curScore;

  writeHierarch(&cur, hc);
  return iters;
}
