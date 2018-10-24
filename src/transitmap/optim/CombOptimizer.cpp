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
int CombOptimizer::optimize(TransitGraph* tg) const {
  LOG(DEBUG) << " == CombOptimizer ==";
  // create optim graph
  OptGraph g(tg, _scorer);

  size_t maxC = maxCard(*g.getNds());
  double solSp = solutionSpaceSize(*g.getNds());
  LOG(DEBUG) << "Optimizing graph of size " << tg->getNodes()->size() <<
    " with max cardinality = " << maxC << " and solution space size = "
    << solSp;

  if (_cfg->createCoreOptimGraph) {
    // TODO: do this exactly as often as M - this should be enough (prove this!)
    for (size_t i = 0; i < 10; i++) {
      g.untangle();
      g.simplify();
      g.split();
    }
  }

  if (_cfg->outOptGraph) {
    util::geo::output::GeoGraphJsonOutput out;
    std::ofstream fstr(_cfg->dbgPath + "/optgraph.json");
    out.print(g, fstr);
  }

  if (_cfg->outputStats) {
    LOG(INFO) << "(stats) Stats for optim graph of '" << tg->getName()
              << std::endl;
    LOG(INFO) << "(stats)   Total node count: " << g.getNumNodes() << " ("
              << g.getNumNodes(true) << " topo, " << g.getNumNodes(false)
              << " non-topo)" << std::endl;
    LOG(INFO) << "(stats)   Total edge count: " << g.getNumEdges() << std::endl;
    LOG(INFO) << "(stats)   Total unique route count: " << g.getNumRoutes()
              << std::endl;
    LOG(INFO) << "(stats)   Max edge route cardinality: "
              << g.getMaxCardinality() << std::endl;
  }

  OrderingConfig c;
  HierarchOrderingConfig hc;


  // iterate over components and optimize all of them separately
  for (const auto nds : util::graph::Algorithm::connectedComponents(g)) {
    optimize(nds, &hc);
  }

  hc.writeFlatCfg(&c);

  Optimizer::expandRelatives(tg, &c);
  tg->setConfig(c);

  return 0;
}

// _____________________________________________________________________________
int CombOptimizer::optimize(const std::set<OptNode*>& g,
                           HierarchOrderingConfig* hc) const {
  size_t maxC = maxCard(g);
  double solSp = solutionSpaceSize(g);
  LOG(DEBUG) << "Optimizing subgraph component of size " << g.size() <<
    " with max cardinality = " << maxC << " and solution space size = "
    << solSp;

  if (maxC == 1) {
    _nullOpt.optimize(g, hc);
  } else {
    _ilpOpt.optimize(g, hc);
  }

  return 0;
}

// _____________________________________________________________________________
size_t CombOptimizer::maxCard(const std::set<OptNode*>& g) {
  size_t ret = 0;
  for (const auto* n : g) {
    for (const auto* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (e->pl().getCardinality() > ret) ret = e->pl().getCardinality();
    }
  }

  return ret;
}

// _____________________________________________________________________________
double CombOptimizer::solutionSpaceSize(const std::set<OptNode*>& g) {
  double ret = 1;
  for (const auto* n : g) {
    for (const auto* e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      ret *= factorial(e->pl().getCardinality());
    }
  }
  return ret;
}
