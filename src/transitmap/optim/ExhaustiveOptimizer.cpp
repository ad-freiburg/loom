// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/optim/ExhaustiveOptimizer.h"

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

  OrderingConfig* c;

  for (OptNode* n : g) {
    for (OptEdge* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;

      // double score = _scorer.getScore( 
    }
  }
  return 0;
}

