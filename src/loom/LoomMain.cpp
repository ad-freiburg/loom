// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "loom/config/ConfigReader.cpp"
#include "loom/config/LoomConfig.h"
#include "loom/optim/CombOptimizer.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "shared/rendergraph/Penalties.h"
#include "shared/rendergraph/RenderGraph.h"
#include "util/geo/PolyLine.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using namespace loom;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  LOGTO(DEBUG, std::cerr) << "Reading graph...";
  shared::rendergraph::RenderGraph g(5, 5);

  g.readFromJson(&std::cin, 3);

  LOGTO(DEBUG, std::cerr) << "Optimizing...";

  double maxCrossPen =
      g.maxDeg() * std::max(cfg.crossPenMultiSameSeg,
                            std::max(cfg.crossPenMultiDiffSeg,
                                     std::max(cfg.stationCrossWeightSameSeg,
                                              cfg.stationCrossWeightDiffSeg)));
  double maxSepPen =
      g.maxDeg() *
      std::max(cfg.separationPenWeight, cfg.stationSeparationWeight);

  // TODO move this into configuration, at least partially
  shared::rendergraph::Penalties pens{maxCrossPen,
                                      maxSepPen,
                                      cfg.crossPenMultiSameSeg,
                                      cfg.crossPenMultiDiffSeg,
                                      cfg.separationPenWeight,
                                      cfg.stationCrossWeightSameSeg,
                                      cfg.stationCrossWeightDiffSeg,
                                      cfg.stationSeparationWeight,
                                      true,
                                      true};

  if (cfg.outputStats) {
    LOGTO(INFO, std::cerr) << "(stats) Stats for graph";
    LOGTO(INFO, std::cerr) << "(stats)   Total node count: "
                           << g.getNds()->size() << " (" << g.numNds(true)
                           << " topo, " << g.numNds(false) << " non-topo)";
    LOGTO(INFO, std::cerr) << "(stats)   Total edge count: " << g.numEdgs();
    LOGTO(INFO, std::cerr) << "(stats)   Total unique line count: "
                           << g.numLines();
    LOGTO(INFO, std::cerr) << "(stats)   Max edge line cardinality: "
                           << g.getMaxLineNum();
    LOGTO(INFO, std::cerr) << "(stats)   Number of poss. solutions: "
                           << g.searchSpaceSize();
    LOGTO(INFO, std::cerr) << "(stats)   Highest node degree: " << g.maxDeg();
  }

  if (cfg.optimMethod == "ilp_impr") {
    optim::ILPEdgeOrderOptimizer ilpEoOptim(&cfg, pens);
    ilpEoOptim.optimize(&g);
  } else if (cfg.optimMethod == "ilp") {
    optim::ILPOptimizer ilpOptim(&cfg, pens);
    ilpOptim.optimize(&g);
  } else if (cfg.optimMethod == "comb") {
    optim::CombOptimizer ilpCombiOptim(&cfg, pens);
    ilpCombiOptim.optimize(&g);
  } else if (cfg.optimMethod == "exhaust") {
    optim::ExhaustiveOptimizer exhausOptim(&cfg, pens);
    exhausOptim.optimize(&g);
  } else if (cfg.optimMethod == "hillc") {
    optim::HillClimbOptimizer hillcOptim(&cfg, pens);
    hillcOptim.optimize(&g);
  } else if (cfg.optimMethod == "anneal") {
    optim::SimulatedAnnealingOptimizer annealOptim(&cfg, pens);
    annealOptim.optimize(&g);
  } else if (cfg.optimMethod == "null") {
    optim::NullOptimizer nullOptim(&cfg, pens);
    nullOptim.optimize(&g);
  }

  util::geo::output::GeoGraphJsonOutput out;
  out.print(g, std::cout);

  return (0);
}
