// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/NullOptimizer.h"
#include "util/Misc.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::NullOptimizer;

// _____________________________________________________________________________
int NullOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                HierarchOrderingConfig* hc,
                                size_t depth) const {
  LOG(DEBUG) << prefix(depth) << "(NullOptimizer) Optimizing component with "
             << g.size() << " nodes.";

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      for (auto etgp : e->pl().etgs) {
        if (etgp.wasCut) continue;
        for (auto ro : e->pl().getRoutes()) {

          for (auto rel : ro.relatives) {
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
  return 0;
}
