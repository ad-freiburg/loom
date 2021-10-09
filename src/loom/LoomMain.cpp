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
#include "loom/optim/GreedyOptimizer.h"
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
  double maxSepPen = g.maxDeg() * std::max(cfg.separationPenWeight,
                                           cfg.stationSeparationWeight);

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

  loom::optim::OptResStats stats;

  if (cfg.optimMethod == "ilp_impr") {
    optim::ILPEdgeOrderOptimizer ilpEoOptim(&cfg, pens);
    stats = ilpEoOptim.optimize(&g);
  } else if (cfg.optimMethod == "ilp") {
    optim::ILPOptimizer ilpOptim(&cfg, pens);
    stats = ilpOptim.optimize(&g);
  } else if (cfg.optimMethod == "comb") {
    optim::CombOptimizer ilpCombiOptim(&cfg, pens);
    stats = ilpCombiOptim.optimize(&g);
  } else if (cfg.optimMethod == "exhaust") {
    optim::ExhaustiveOptimizer exhausOptim(&cfg, pens);
    stats = exhausOptim.optimize(&g);
  } else if (cfg.optimMethod == "hillc") {
    optim::HillClimbOptimizer hillcOptim(&cfg, pens, false);
    stats = hillcOptim.optimize(&g);
  } else if (cfg.optimMethod == "hillc-random") {
    optim::HillClimbOptimizer hillcOptim(&cfg, pens, true);
    stats = hillcOptim.optimize(&g);
  } else if (cfg.optimMethod == "anneal") {
    optim::SimulatedAnnealingOptimizer annealOptim(&cfg, pens, false);
    stats = annealOptim.optimize(&g);
  } else if (cfg.optimMethod == "anneal-random") {
    optim::SimulatedAnnealingOptimizer annealOptim(&cfg, pens, true);
    stats = annealOptim.optimize(&g);
  } else if (cfg.optimMethod == "greedy") {
    optim::GreedyOptimizer greedyOptim(&cfg, pens, false);
    stats = greedyOptim.optimize(&g);
  } else if (cfg.optimMethod == "greedy-lookahead") {
    optim::GreedyOptimizer greedyOptim(&cfg, pens, true);
    stats = greedyOptim.optimize(&g);
  } else if (cfg.optimMethod == "null") {
    optim::NullOptimizer nullOptim(&cfg, pens);
    stats = nullOptim.optimize(&g);
  } else {
    LOG(ERROR) << "Unknown optimization method " << cfg.optimMethod
               << std::endl;
    exit(1);
  }

  util::geo::output::GeoGraphJsonOutput out;

  if (cfg.outputStats) {
    util::json::Dict jsonStats = {
        {"statistics",
         util::json::Dict{
             {"input_num_nodes", stats.numNodes},
             {"input_num_edges", stats.numEdges},
             {"input_max_number_lines", stats.maxLineCard},
             {"input_solution_space_size", stats.solutionSpaceSize},
             {"input_nontrivial_comps", stats.nonTrivialComponents},
             {"input_nontrivial_comps_searchspace_one",
              stats.numCompsSolSpaceOne},
             {"input_max_num_nodes_in_comps", stats.maxNumNodesPerComp},
             {"input_max_num_edges_in_comps", stats.maxNumEdgesPerComp},
             {"input_max_number_lines_in_comps", stats.maxCardPerComp},
             {"input_max_solution_space_size_in_comps", stats.maxCompSolSpace},
             {"runs", stats.runs},
             {"avg_solve_time", stats.avgSolveTime},
             {"avg_iterations", stats.avgIterations},
             {"avg_score", stats.avgScore},
             {"avg_num_same_seg_crossings", stats.avgSameSegCross},
             {"avg_num_diff_seg_crossings", stats.avgDiffSegCross},
             {"avg_num_separations", stats.avgSeps},
             {"best_num_same_seg_crossings", stats.sameSegCrossings},
             {"best_num_diff_seg_crossings", stats.diffSegCrossings},
             {"best_num_separations", stats.separations},
             {"best_score", stats.score}}}};
    out.print(g, std::cout, jsonStats);
  } else {
    out.print(g, std::cout);
  }

  return (0);
}
