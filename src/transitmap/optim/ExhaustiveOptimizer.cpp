// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include "transitmap/optim/ExhaustiveOptimizer.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::ExhaustiveOptimizer;

// _____________________________________________________________________________
int ExhaustiveOptimizer::optimize(TransitGraph* tg) const {
  // TODO

  return 0;
}

// _____________________________________________________________________________
int ExhaustiveOptimizer::optimize(const std::set<OptNode*>& g,
                                  HierarchOrderingConfig* hc) const {
  OptOrderingConfig best;
  OptOrderingConfig cur;
  double bestScore = DBL_MAX;

  // fixed order list of optim graph edges
  std::vector<OptEdge*> edges;

  for (auto n : g) {
    for (auto e : n->getAdjList()) {
      if (n == e->getFrom()) edges.push_back(e);
    }
  }

  // this guarantees that all the orderings are sorted!
  initialConfig(g, &cur);

  size_t iters = 0;
  bool running = true;

  while (running) {
    double score = _optScorer.getScore(g, cur);
    // LOG(DEBUG) << "Iteration " << iters << ", current score: " << score;
    iters++;

    if (score < bestScore) {
      bestScore = score;
      best = cur;
    }

    if (score == 0) {
      LOG(DEBUG) << "Found optimal score 0 prematurely after " << iters << " iterations!";
      writeHierarch(&best, hc);
      return 0;
    }

    for (size_t i = 0; i < edges.size(); i++) {
      if (std::next_permutation(cur[edges[i]].begin(), cur[edges[i]].end())) {
        break;
      } else if (i == edges.size() - 1) {
        running = false;
      } else {
        // reset
        std::sort(cur[edges[i]].begin(), cur[edges[i]].end());
      }
    }

  }

  LOG(DEBUG) << "Found optimal score " << bestScore << " after " << iters << " iterations!";

  writeHierarch(&best, hc);
  return 0;
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderingConfig* cfg) const {
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      (*cfg)[e] = std::vector<const graph::Route*>(e->pl().getCardinality());
      for (size_t i = 0; i < (*cfg)[e].size(); i++)
        (*cfg)[e][i] = e->pl().getRoutes()[i].route;

      std::sort((*cfg)[e].begin(), (*cfg)[e].end());
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

          if (ro.route->relativeTo()) continue;

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
