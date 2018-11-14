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
#include "util/graph/Algorithm.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace optim;
using namespace transitmapper::graph;
using transitmapper::optim::CombOptimizer;

// _____________________________________________________________________________
int CombOptimizer::optimizeComp(const std::set<OptNode*>& g,
                           HierarchOrderingConfig* hc) const {
  size_t maxC = maxCard(g);
  double solSp = solutionSpaceSize(g);

  if (maxC == 1) {
    LOG(DEBUG) << "(Null optimizer)";
    _nullOpt.optimizeComp(g, hc);
  } else if (solSp < 5000) {
    LOG(DEBUG) << "(Exhaustive optimizer)";
    _exhausOpt.optimizeComp(g, hc);
  } else {
    LOG(DEBUG) << "(ILP optimizer)";
    _ilpOpt.optimizeComp(g, hc);
  }

  return 0;
}

