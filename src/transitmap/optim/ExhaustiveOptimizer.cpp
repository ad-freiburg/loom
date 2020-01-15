// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "shared/linegraph/Route.h"
#include "transitmap/optim/ExhaustiveOptimizer.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::ExhaustiveOptimizer;
using shared::linegraph::Route;

// _____________________________________________________________________________
int ExhaustiveOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                      HierarOrderCfg* hc,
                                      size_t depth) const {
  LOG(DEBUG) << prefix(depth)
             << "(ExhaustiveOptimizer) Optimizing component with " << g.size()
             << " nodes.";

  OptOrderCfg best, cur, null;
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

  double curScore = _optScorer.getCrossingScore(og, g, cur);
  if (_cfg->splittingOpt) curScore += _optScorer.getSplittingScore(og, g, cur);

  bestScore = curScore;
  best = cur;

  while (true) {
    if (bestScore == 0) {
      LOG(DEBUG) << prefix(depth) << "Found optimal score 0 prematurely after "
                 << iters << " iterations!";
      writeHierarch(&best, hc);
      return 0;
    }

    iters++;

    if (iters - last == 10000) {
      LOG(DEBUG) << prefix(depth) << "@ " << iters;
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

    double curScore = _optScorer.getCrossingScore(og, g, cur);
    if (_cfg->splittingOpt) curScore += _optScorer.getSplittingScore(og, g, cur);

    if (curScore < bestScore) {
      bestScore = curScore;
      best = cur;
    }
  }

  LOG(DEBUG) << prefix(depth) << "Found optimal score " << bestScore
             << " after " << iters << " iterations!";

  writeHierarch(&best, hc);

  return iters;
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderCfg* cfg) const {
  initialConfig(g, cfg, false);
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderCfg* cfg,
                                        bool sorted) const {
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      (*cfg)[e] = std::vector<const Route*>(e->pl().getCardinality());
      size_t p = 0;
      for (size_t i = 0; i < e->pl().getRoutes().size(); i++) {
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
void ExhaustiveOptimizer::writeHierarch(OptOrderCfg* cfg,
                                        HierarOrderCfg* hc) const {
  for (auto ep : *cfg) {
    auto e = ep.first;

    for (auto etgp : e->pl().etgs) {
      if (etgp.wasCut) continue;
      for (auto r : ep.second) {
        // get the corresponding route occurance in the opt graph edge
        // TODO: replace this as soon as a lookup function is present in OptRO
        OptRO optRO;
        for (auto ro : e->pl().getRoutes()) {
          if (r == ro.route) optRO = ro;
        }

        for (auto rel : optRO.relatives) {
          // retrieve the original route pos
          size_t p = etgp.etg->pl().getRoutePos(rel);
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
