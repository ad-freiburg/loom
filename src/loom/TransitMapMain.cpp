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
#include "loom/config/TransitMapConfig.h"
#include "loom/graph/GraphBuilder.h"
#include "loom/graph/Penalties.h"
#include "loom/graph/RenderGraph.h"
#include "loom/optim/CombOptimizer.h"
#include "loom/optim/ILPEdgeOrderOptimizer.h"
#include "loom/optim/Scorer.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace loom;
using namespace util::geo;
using std::string;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  LOG(INFO) << "Reading graph...";
  loom::graph::RenderGraph g(cfg.lineWidth, cfg.lineSpacing);
  loom::graph::GraphBuilder b(&cfg);

  g.readFromJson(&std::cin, cfg.inputSmoothing);

  LOG(INFO) << "Creating node fronts...";
  b.writeMainDirs(&g);

  LOG(INFO) << "Optimizing...";

  double maxCrossPen =
      g.maxDeg() * (cfg.crossPenMultiSameSeg > cfg.crossPenMultiDiffSeg
                              ? cfg.crossPenMultiSameSeg
                              : cfg.crossPenMultiDiffSeg);
  double maxSplitPen = g.maxDeg() * cfg.splitPenWeight;

  // TODO move this into configuration, at least partially
  loom::graph::Penalties pens{maxCrossPen,
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
    LOG(INFO) << "(stats) Stats for graph";
    LOG(INFO) << "(stats)   Total node count: " << g.getNds()->size() << " ("
              << g.getNumNds(true) << " topo, " << g.getNumNds(false)
              << " non-topo)";
    LOG(INFO) << "(stats)   Total edge count: " << g.numEdgs();
    LOG(INFO) << "(stats)   Total unique line count: " << g.getNumLines();
    LOG(INFO) << "(stats)   Max edge line cardinality: "
              << g.getMaxLineNum();
    LOG(INFO) << "(stats)   Number of poss. solutions: "
              << scorer.getNumPossSolutions();
    LOG(INFO) << "(stats)   Highest node degree: " << g.maxDeg();
  }

  LOG(INFO) << "(stats) Max crossing pen: " << maxCrossPen;
  LOG(INFO) << "(stats) Max splitting pen: " << maxSplitPen;

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

  LOG(INFO) << "(stats) Total graph score AFTER optim is -- "
            << scorer.getScore() << " -- (excl. unavoidable crossings!)";
  LOG(INFO) << "(stats)   Per node graph score: "
            << scorer.getScore() / g.getNds()->size();
  LOG(INFO) << "(stats)   Crossings: " << scorer.getNumCrossings()
            << " (score: " << scorer.getCrossScore() << ")";
  LOG(INFO) << "(stats)   Separations: " << scorer.getNumSeparations()
            << " (score: " << scorer.getSeparationScore() << ")";

  return (0);
}
