// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <unordered_map>
#include "loom/optim/OracleOptimizer.h"
#include "shared/linegraph/Line.h"
#include "util/log/Log.h"

using namespace loom;
using namespace optim;
using loom::optim::OracleOptimizer;
using shared::linegraph::Line;
using shared::rendergraph::HierarOrderCfg;

// _____________________________________________________________________________
int OracleOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                  HierarOrderCfg* hc, size_t depth) const {
  UNUSED(depth);
  OptOrderCfg cur, null;

  // fixed order list of optim graph edges
  std::vector<OptEdge*> edges;

  for (auto n : g)
    for (auto e : n->getAdjList())
      if (n == e->getFrom() && e->pl().getCardinality() > 1) edges.push_back(e);

  // this guarantees that all the orderings are sorted!
  initialConfig(g, &null, true);

  writeHierarch(&cur, hc);
  return 0;
}
