// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "transitmap/config/ConfigReader.cpp"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/GraphBuilder.h"
#include "transitmap/graph/Node.h"
#include "transitmap/graph/Penalties.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/optim/CombOptimizer.h"
#include "transitmap/optim/ILPEdgeOrderOptimizer.h"
#include "transitmap/optim/Scorer.h"
#include "transitmap/output/SvgOutput.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace transitmapper;
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
  transitmapper::graph::TransitGraph g;
  transitmapper::graph::GraphBuilder b(&cfg);

  if (!b.build(&(std::cin), &g)) {
    exit(1);
  }

  LOG(INFO) << "Creating node fronts...";
  b.writeMainDirs(&g);

  if (cfg.dontExpandStations) {
    b.writeStationGeoms(&g);
  }

  if (cfg.expandFronts) {
    b.expandOverlappinFronts(&g);
  }

  if (!cfg.dontExpandStations) {
    b.writeStationGeoms(&g);
  }

  b.createMetaNodes(&g);

  LOG(INFO) << "Optimizing...";

  double maxCrossPen =
      g.getMaxDegree() * (cfg.crossPenMultiSameSeg > cfg.crossPenMultiDiffSeg
                              ? cfg.crossPenMultiSameSeg
                              : cfg.crossPenMultiDiffSeg);
  double maxSplitPen = g.getMaxDegree() * cfg.splitPenWeight;

  // TODO move this into configuration, at least partially
  transitmapper::graph::Penalties pens{maxCrossPen,
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
              << g.getNumNodes(true) << " topo, " << g.getNumNodes(false)
              << " non-topo)";
    LOG(INFO) << "(stats)   Total edge count: " << g.getNumEdges();
    LOG(INFO) << "(stats)   Total unique route count: " << g.getNumRoutes();
    LOG(INFO) << "(stats)   Max edge route cardinality: "
              << g.getMaxCardinality();
    LOG(INFO) << "(stats)   Number of poss. solutions: "
              << scorer.getNumPossSolutions();
    LOG(INFO) << "(stats)   Highest node degree: " << g.getMaxDegree();
  }

  LOG(INFO) << "(stats) Max crossing pen: " << maxCrossPen;
  LOG(INFO) << "(stats) Max splitting pen: " << maxSplitPen;

  if (!cfg.noOptim) {
    if (cfg.optimMethod == "ilp_impr") {
      LOG(INFO) << "(ILP impr optimizer)";
      optim::ILPEdgeOrderOptimizer ilpEoOptim(&cfg, &scorer);
      ilpEoOptim.optimize(&g);
    } else if (cfg.optimMethod == "ilp") {
      LOG(INFO) << "(ILP optimizer)";
      optim::ILPOptimizer ilpOptim(&cfg, &scorer);
      ilpOptim.optimize(&g);
    } else if (cfg.optimMethod == "comb") {
      LOG(INFO) << "(Comb optimizer)";
      optim::CombOptimizer ilpCombiOptim(&cfg, &scorer);
      ilpCombiOptim.optimize(&g);
    } else if (cfg.optimMethod == "exhaust") {
      LOG(INFO) << "(Exhaustive optimizer)";
      optim::ExhaustiveOptimizer exhausOptim(&cfg, &scorer);
      exhausOptim.optimize(&g);
    } else if (cfg.optimMethod == "hillc") {
      LOG(INFO) << "(Hillclimbing optimizer)";
      optim::HillClimbOptimizer hillcOptim(&cfg, &scorer);
      hillcOptim.optimize(&g);
    } else if (cfg.optimMethod == "anneal") {
      LOG(INFO) << "(Simulated annealing optimizer)";
      optim::SimulatedAnnealingOptimizer annealOptim(&cfg, &scorer);
      annealOptim.optimize(&g);
    }

    LOG(INFO) << "(stats) Total graph score AFTER optim is -- "
              << scorer.getScore() << " -- (excl. unavoidable crossings!)";
    LOG(INFO) << "(stats)   Per node graph score: "
              << scorer.getScore() / g.getNds()->size();
    LOG(INFO) << "(stats)   Crossings: " << scorer.getNumCrossings()
              << " (score: " << scorer.getCrossScore() << ")";
    LOG(INFO) << "(stats)   Separations: " << scorer.getNumSeparations()
              << " (score: " << scorer.getSeparationScore() << ")";
  }

  if (cfg.renderMethod == "svg") {
    std::string path = cfg.outputPath;
    LOG(INFO) << "Outputting to SVG " << path << " ..." << std::endl;
    std::ofstream o;
    o.open(path);
    output::SvgOutput svgOut(&o, &cfg, &scorer);
    svgOut.print(g);
  }

  if (cfg.renderMethod == "svg_sep") {
    {
      std::string path = cfg.outputPath + "/edges.svg";
      cfg.renderEdges = 1;
      cfg.renderStations = 0;
      cfg.renderNodeConnections = 0;
      LOG(INFO) << "Outputting edge SVG to " << path << " ..." << std::endl;
      std::ofstream o;
      o.open(path);
      output::SvgOutput svgOut(&o, &cfg, &scorer);
      svgOut.print(g);
    }

    {
      std::string path = cfg.outputPath + "/nodes.svg";
      cfg.renderEdges = 0;
      cfg.renderStations = 0;
      cfg.renderNodeConnections = 1;
      LOG(INFO) << "Outputting node SVG to " << path << " ..." << std::endl;
      std::ofstream o;
      o.open(path);
      output::SvgOutput svgOut(&o, &cfg, &scorer);
      svgOut.print(g);
    }

    {
      std::string path = cfg.outputPath + "/stations.svg";
      cfg.renderEdges = 0;
      cfg.renderStations = 1;
      cfg.renderNodeConnections = 0;
      LOG(INFO) << "Outputting station SVG to " << path << " ..." << std::endl;
      std::ofstream o;
      o.open(path);
      output::SvgOutput svgOut(&o, &cfg, &scorer);
      svgOut.print(g);
    }
  }

  if (!cfg.worldFilePath.empty()) {
    LOG(INFO) << "Writing world file for this map to " << cfg.worldFilePath;

    std::ofstream file;
    file.open(cfg.worldFilePath);
    if (file) {
      file << 1 / cfg.outputResolution << std::endl
           << 0 << std::endl
           << 0 << std::endl
           << -1 / cfg.outputResolution << std::endl
           << std::fixed
           << g.getBoundingBox(cfg.outputPadding).getLowerLeft().getX()
           << std::endl
           << g.getBoundingBox(cfg.outputPadding).getUpperRight().getY()
           << std::endl;
      file.close();
    }
  }
  return (0);
}
