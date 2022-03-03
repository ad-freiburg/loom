// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "loom/optim/ExhaustiveOptimizer.h"
#include "shared/linegraph/Line.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using loom::optim::ExhaustiveOptimizer;
using shared::linegraph::Line;
using shared::rendergraph::HierarOrderCfg;

// _____________________________________________________________________________
double ExhaustiveOptimizer::optimizeComp(OptGraph* og,
                                         const std::set<OptNode*>& g,
                                         HierarOrderCfg* hc, size_t depth,
                                         OptResStats& stats) const {
  UNUSED(og);
  UNUSED(stats);
  LOGTO(DEBUG, std::cerr) << prefix(depth)
                          << "(ExhaustiveOptimizer) Optimizing component with "
                          << g.size() << " nodes.";

  T_START(1);

  OptOrderCfg best, cur, null;
  double bestScore = DBL_MAX;

  // fixed order list of optim graph edges
  std::vector<OptEdge*> edges;

  for (auto n : g)
    for (auto e : n->getAdjList())
      if (n == e->getFrom()) edges.push_back(e);

  // this guarantees that all the orderings are sorted, which we need for
  // std::next_permutation below!
  initialConfig(g, &null, true);
  cur = null;

  double iters = 0;
  double last = 0;
  bool running = true;

  double curScore = _optScorer.getCrossingScore(g, cur);
  if (_optScorer.optimizeSep())
    curScore += _optScorer.getSeparationScore(g, cur);

  bestScore = curScore;
  best = cur;

  double solSp = solutionSpaceSize(g);

  // don't try if it is pointless, assuming we can make 100.000
  // iterations per second
  if ((solSp / 50000) > (60 * 60 * 6)) {
    std::stringstream ss;
    ss << "Exhaustive search would take too long (over "
       << ((solSp / 50000) / (60 * 60))
       << " hours even assuming we can check 50.000 configurations per second";
    throw std::runtime_error(ss.str());
  }

  double itTime = 0;

  while (true) {
    T_START(iter);
    if (bestScore == 0) {
      LOGTO(DEBUG, std::cerr)
          << prefix(depth) << "Found optimal score 0 prematurely after "
          << iters << " iterations!";
      writeHierarch(&best, hc);
      return 0;
    }

    iters++;

    if (fabs((iters - last) - 10000) < 1) {
      LOGTO(DEBUG, std::cerr)
          << prefix(depth) << "@ " << iters << "/" << solSp << " ("
          << int(((1.0 * iters) / (1.0 * solSp)) * 100) << "%, "
          << (10000 / (itTime / 1000)) << " its/s)";
      last = iters;
      itTime = 0;
    }

    for (size_t i = 0; i < edges.size(); i++) {
      if (std::next_permutation(cur[edges[i]].begin(), cur[edges[i]].end())) {
        break;
      } else if (i == edges.size() - 1) {
        running = false;
      }
    }

    if (!running) break;

    if (_optScorer.optimizeSep())
      curScore = _optScorer.getTotalScore(g, cur);
    else
      curScore = _optScorer.getCrossingScore(g, cur);

    if (curScore < bestScore) {
      bestScore = curScore;
      best = cur;
    }
    itTime += T_STOP(iter);
  }

  LOGTO(DEBUG, std::cerr) << prefix(depth) << "Found optimal score "
                          << bestScore << " after " << iters << " iterations!";

  writeHierarch(&best, hc);

  return T_STOP(1);
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderCfg* cfg) const {
  initialConfig(g, cfg, false);
}

// _____________________________________________________________________________
void ExhaustiveOptimizer::initialConfig(const std::set<OptNode*>& g,
                                        OptOrderCfg* cfg, bool sorted) const {
  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      (*cfg)[e] = std::vector<const Line*>(e->pl().getCardinality());
      size_t p = 0;
      for (size_t i = 0; i < e->pl().getLines().size(); i++) {
        (*cfg)[e][p] = e->pl().getLines()[i].line;
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

    for (auto lnEdgPart : e->pl().lnEdgParts) {
      if (lnEdgPart.wasCut) continue;
      for (auto r : ep.second) {
        // get the corresponding route occurance in the opt graph edge
        // TODO: replace this as soon as a lookup function is present in OptLO
        OptLO optRO;
        for (auto ro : e->pl().getLines()) {
          if (r == ro.line) optRO = ro;
        }

        for (auto rel : optRO.relatives) {
          // retrieve the original line pos
          size_t p = lnEdgPart.lnEdg->pl().linePos(rel);
          if (!(lnEdgPart.dir ^ e->pl().lnEdgParts.front().dir)) {
            (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].insert(
                (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].begin(), p);
          } else {
            (*hc)[lnEdgPart.lnEdg][lnEdgPart.order].push_back(p);
          }
        }
      }
    }
  }
}
