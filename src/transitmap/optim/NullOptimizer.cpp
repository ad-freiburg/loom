// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/NullOptimizer.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::NullOptimizer;

// _____________________________________________________________________________
int NullOptimizer::optimizeComp(const std::set<OptNode*>& g,
                           HierarchOrderingConfig* hc) const {

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      assert(e->pl().getCardinality() == 1);
      for (auto etgp : e->pl().etgs) {
        if (etgp.wasCut) continue;
        for (auto ro : e->pl().getRoutes()) {
          assert((*hc)[etgp.etg][etgp.order].size() == 0);

          for (auto rel : ro.relatives) {
            // retrieve the original route pos
            size_t p = etgp.etg->getRoutePos(rel);
            (*hc)[etgp.etg][etgp.order].push_back(p);
          }
        }
      }
    }
  }
  return 0;
}
