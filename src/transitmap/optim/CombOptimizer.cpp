// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <glpk.h>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <thread>
#include "transitmap/graph/OrderingConfig.h"
#include "transitmap/optim/CombOptimizer.h"
#include "transitmap/optim/OptGraph.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Algorithm.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::CombOptimizer;

// _____________________________________________________________________________
int CombOptimizer::optimizeComp(OptGraph* og, const std::set<OptNode*>& g,
                                HierarchOrderingConfig* hc, size_t depth) const {
  size_t maxC = maxCard(g);
  double solSp = solutionSpaceSize(g);

  LOG(DEBUG) << prefix(depth) << "(CombOptimizer) Optimizing comp with " << g.size()
             << " nodes, max card " << maxC << ", sol space size "
             << solSp;

  if (maxC == 1) {
    _nullOpt.optimizeComp(og, g, hc, depth + 1);
  } else if (solSp < 10) {
    _exhausOpt.optimizeComp(og, g, hc, depth + 1);
  } else {
    _ilpOpt.optimizeComp(og, g, hc, depth + 1);
  }

  return 0;
}
