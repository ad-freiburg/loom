// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "transitmap/optim/ExhaustiveOptimizer.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::ExhaustiveOptimizer;

// _____________________________________________________________________________
int ExhaustiveOptimizer::optimizeComp(const std::set<OptNode*>& g,
                                      HierarchOrderingConfig* hc) const {
  OptOrderingConfig best, cur, null;
  double bestScore = DBL_MAX;

  // fixed order list of optim graph edges
  std::vector<OptEdge*> edges;

  for (auto n : g)
    for (auto e : n->getAdjList())
      if (n == e->getFrom()) edges.push_back(e);

  // this guarantees that all the orderings are sorted!
  initialConfig(g, &null, true);
  cur = null;

  size_t iters = 0;
  size_t last = 0;
  bool running = true;

  double curScore = _optScorer.getCrossingScore(g, cur);
  if (_cfg->splittingOpt) curScore += _optScorer.getSplittingScore(g, cur);

  bestScore = curScore;
  best = cur;

  while (true) {
    if (bestScore == 0) {
      LOG(DEBUG) << "Found optimal score 0 prematurely after " << iters
                 << " iterations!";
      writeHierarch(&best, hc);
      return 0;
    }

    iters++;

    if (iters - last == 10000) {
      LOG(DEBUG) << "@ " << iters;
      last = iters;
    }

    for (size_t i = 0; i < edges.size(); i++) {
      if (std::next_permutation(cur[edges[i]].begin(), cur[edges[i]].end())) {
        break;
      } else if (i == edges.size() - 1) {
        running = false;
      } else {
        // reset
        cur[edges[i]] = null[edges[i]];
      }
    }

    if (!running) break;

    double curScore = _optScorer.getCrossingScore(g, cur);
    if (_cfg->splittingOpt) curScore += _optScorer.getSplittingScore(g, cur);

    if (curScore < bestScore) {
      bestScore = curScore;
      best = cur;
    }
  }

  LOG(DEBUG) << "Found optimal score " << bestScore << " after " << iters
             << " iterations!";

  writeHierarch(&best, hc);
  return iters;
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderingConfig* cfg) const {
  initialConfig(g, cfg, false);
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderingConfig* cfg,
                                        bool sorted) const {
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      (*cfg)[e] = std::vector<const graph::Route*>(e->pl().getCardinality());
      size_t p = 0;
      for (size_t i = 0; i < e->pl().getRoutes().size(); i++) {
        if (e->pl().getRoutes()[i].relativeTo) continue;
        (*cfg)[e][p] = e->pl().getRoutes()[i].route;
        p++;
      }

      if (sorted) {
        std::sort((*cfg)[e].begin(), (*cfg)[e].end());
      } else {
        std::random_shuffle((*cfg)[e].begin(), (*cfg)[e].end());
      }
    }
  }
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::writeHierarch(OptOrderingConfig* cfg,
                                        HierarchOrderingConfig* hc) const {
  for (auto ep : *cfg) {
    auto e = ep.first;

    for (auto etgp : e->pl().etgs) {
      if (etgp.wasCut) continue;
      for (auto r : ep.second) {
        for (size_t p = 0; p < etgp.etg->getCardinality(); p++) {
          auto ro = (*etgp.etg->getRoutes())[p];
          if (!(r == ro.route)) continue;

          if (std::find(e->pl().getRoutes().begin(), e->pl().getRoutes().end(),
                        ro) == e->pl().getRoutes().end())
            continue;

          // if (ro.route->relativeTo()) continue;

          if (!(etgp.dir ^ e->pl().etgs.front().dir)) {
            (*hc)[etgp.etg][etgp.order].insert(
                (*hc)[etgp.etg][etgp.order].begin(), p);
          } else {
            (*hc)[etgp.etg][etgp.order].push_back(p);
          }
        }
      }
    }
  }
}
