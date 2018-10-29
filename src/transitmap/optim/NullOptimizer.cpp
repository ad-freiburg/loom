// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/NullOptimizer.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::NullOptimizer;

// _____________________________________________________________________________
int NullOptimizer::optimize(TransitGraph* tg) const {
  // TODO

  return 0;
}

// _____________________________________________________________________________
int NullOptimizer::optimize(const std::set<OptNode*>& g,
                           HierarchOrderingConfig* hc) const {

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      assert(e->pl().getCardinality() == 1);
      for (auto etgp : e->pl().etgs) {
        if (etgp.wasCut) continue;
        for (size_t p = 0; p < etgp.etg->getCardinality(); p++) {
          auto ro = (*etgp.etg->getRoutes())[p];

          if (std::find(e->pl().getRoutes().begin(),
                        e->pl().getRoutes().end(),
                        ro) == e->pl().getRoutes().end())
            continue;

          if (ro.route->relativeTo()) continue;

          assert((*hc)[etgp.etg][etgp.order].size() == 0);
          (*hc)[etgp.etg][etgp.order].push_back(p);
        }
      }
    }
  }
  return 0;
}

