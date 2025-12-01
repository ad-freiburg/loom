// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <chrono>
#include <cstdio>
#include <fstream>
#include <thread>

#include "loom/optim/CombOptimizer.h"
#include "loom/optim/OptGraph.h"
#include "shared/rendergraph/OrderCfg.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Algorithm.h"
#include "util/log/Log.h"

using loom::optim::CombOptimizer;
using shared::rendergraph::HierarOrderCfg;
using util::DEBUG;

// _____________________________________________________________________________
double CombOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                   HierarOrderCfg* hc, size_t depth,
                                   OptResStats& stats) const {
  size_t maxC = maxCard(g);
  double solSp = solutionSpaceSize(g);

  LOGTO(DEBUG, std::cerr) << prefix(depth)
                          << "(CombOptimizer) Optimizing comp with " << g.size()
                          << " nodes, max card " << maxC << ", sol space size "
                          << solSp;

  if (maxC == 1) {
    return _nullOpt.optimizeComp(og, g, hc, depth + 1, stats);
  } else if (solSp < 500) {
    return _exhausOpt.optimizeComp(og, g, hc, depth + 1, stats);
  } else {
    if (_forceILP) return _ilpOpt.optimizeComp(og, g, hc, depth + 1, stats);
#if defined GUROBI_FOUND || defined GLPK_FOUND || defined COIN_FOUND
    return _ilpOpt.optimizeComp(og, g, hc, depth + 1, stats);
#else
    return _hillcOpt.optimizeComp(og, g, hc, depth + 1, stats);
#endif
  }
}
