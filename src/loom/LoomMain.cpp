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
#include "loom/graph/GraphBuilder.h"
#include "loom/optim/CombOptimizer.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/Scorer.h"
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
  loom::graph::GraphBuilder b(&cfg);

  g.readFromJson(&std::cin, 3);

  LOGTO(DEBUG, std::cerr) << "Creating node fronts...";
  b.writeMainDirs(&g);

  LOGTO(DEBUG, std::cerr) << "Optimizing...";

  double maxCrossPen =
      g.maxDeg() * (cfg.crossPenMultiSameSeg > cfg.crossPenMultiDiffSeg
                        ? cfg.crossPenMultiSameSeg
                        : cfg.crossPenMultiDiffSeg);
  double maxSplitPen = g.maxDeg() * cfg.splitPenWeight;

  // TODO move this into configuration, at least partially
  shared::rendergraph::Penalties pens{maxCrossPen,
                                      maxSplitPen,
                                      cfg.crossPenMultiSameSeg,
                                      cfg.crossPenMultiDiffSeg,
                                      cfg.splitPenWeight,
                                      cfg.stationCrossWeightSameSeg,
                                      cfg.stationCrossWeightDiffSeg,
                                      cfg.stationSplitWeight,
                                      true,
                                      true};

  optim::Scorer scorer(&g, pens);

  if (cfg.outputStats) {
    LOGTO(INFO, std::cerr) << "(stats) Stats for graph";
    LOGTO(INFO, std::cerr) << "(stats)   Total node count: "
                           << g.getNds()->size() << " (" << g.getNumNds(true)
                           << " topo, " << g.getNumNds(false) << " non-topo)";
    LOGTO(INFO, std::cerr) << "(stats)   Total edge count: " << g.numEdgs();
    LOGTO(INFO, std::cerr) << "(stats)   Total unique line count: "
                           << g.getNumLines();
    LOGTO(INFO, std::cerr) << "(stats)   Max edge line cardinality: "
                           << g.getMaxLineNum();
    LOGTO(INFO, std::cerr) << "(stats)   Number of poss. solutions: "
                           << scorer.getNumPossSolutions();
    LOGTO(INFO, std::cerr) << "(stats)   Highest node degree: " << g.maxDeg();
  }

  if (cfg.optimMethod == "ilp_impr") {
    optim::ILPEdgeOrderOptimizer ilpEoOptim(&cfg, &scorer);
    ilpEoOptim.optimize(&g);
  } else if (cfg.optimMethod == "ilp") {
    optim::ILPOptimizer ilpOptim(&cfg, &scorer);
    ilpOptim.optimize(&g);
  } else if (cfg.optimMethod == "comb") {
    optim::CombOptimizer ilpCombiOptim(&cfg, &scorer);
    ilpCombiOptim.optimize(&g);
  } else if (cfg.optimMethod == "exhaust") {
    optim::ExhaustiveOptimizer exhausOptim(&cfg, &scorer);
    exhausOptim.optimize(&g);
  } else if (cfg.optimMethod == "hillc") {
    optim::HillClimbOptimizer hillcOptim(&cfg, &scorer);
    hillcOptim.optimize(&g);
  } else if (cfg.optimMethod == "anneal") {
    optim::SimulatedAnnealingOptimizer annealOptim(&cfg, &scorer);
    annealOptim.optimize(&g);
  } else if (cfg.optimMethod == "null") {
    optim::NullOptimizer nullOptim(&cfg, &scorer);
    nullOptim.optimize(&g);
  }

  // LOGTO(INFO,std::cerr) << "(stats) Total graph score AFTER optim is -- "
  // << scorer.getScore() << " -- (excl. unavoidable crossings!)";
  // LOGTO(INFO,std::cerr) << "(stats)   Per node graph score: "
  // << scorer.getScore() / g.getNds()->size();
  // LOGTO(INFO,std::cerr) << "(stats)   Crossings: " <<
  // scorer.getNumCrossings()
  // << " (score: " << scorer.getCrossScore() << ")";
  // LOGTO(INFO,std::cerr) << "(stats)   Separations: " <<
  // scorer.getNumSeparations()
  // << " (score: " << scorer.getSeparationScore() << ")";
  util::geo::output::GeoGraphJsonOutput out;
  out.print(g, std::cout);

  return (0);
}
