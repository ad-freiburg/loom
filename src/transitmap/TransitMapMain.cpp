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

// new
#include "transitmap/graph/RenderLineGraph.h"

using transitmapper::config::Config;
using namespace util::geo;
using std::string;

// new
using transitmapper::graph::RenderLineGraph;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  Config cfg;

  ConfigReader cr;
  cr.read(&cfg, argc, argv);

  // ______________________________________________________ new
  transitmapper::graph::GraphBuilder bn(&cfg);

  RenderLineGraph tg(&cfg);
  tg.readFromJson(&(std::cin));

  // TODO: implement combinePartnerRoutes in optim graph build

  LOG(INFO) << "Creating node fronts...";
  bn.writeFronts(&tg);

  if (cfg.dontExpandStations) {
    tg.writeStatGeoms();
  }

  exit(0);

  // ______________________________________________________ /new

  LOG(INFO) << "Reading graph...";
  transitmapper::graph::TransitGraph g(cfg.name, cfg.projectionString);
  transitmapper::graph::GraphBuilder b(&cfg);

  // read the graph from std::cin to g
  if (!b.build(&(std::cin), &g)) {
    exit(1);
  }

  if (cfg.collapseLinePartners) {
    b.combinePartnerRoutes(&g);
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

  transitmapper::optim::Scorer scorer(&g, pens);

  if (cfg.outputStats) {
    LOG(INFO) << "(stats) Stats for graph '" << g.getName();
    LOG(INFO) << "(stats)   Total node count: " << g.getNumNodes() << " ("
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
      LOG(DEBUG) << "(ILP impr optimizer)";
      transitmapper::optim::ILPEdgeOrderOptimizer ilpEoOptim(&cfg, &scorer);
      ilpEoOptim.optimize(&g);
    } else if (cfg.optimMethod == "ilp") {
      LOG(DEBUG) << "(ILP optimizer)";
      transitmapper::optim::ILPOptimizer ilpOptim(&cfg, &scorer);
      ilpOptim.optimize(&g);
    } else if (cfg.optimMethod == "comb") {
      LOG(DEBUG) << "(Comb optimizer)";
      transitmapper::optim::CombOptimizer ilpCombiOptim(&cfg, &scorer);
      ilpCombiOptim.optimize(&g);
    } else if (cfg.optimMethod == "exhaust") {
      LOG(DEBUG) << "(Exhaustive optimizer)";
      transitmapper::optim::ExhaustiveOptimizer exhausOptim(&cfg, &scorer);
      exhausOptim.optimize(&g);
    } else if (cfg.optimMethod == "hillc") {
      LOG(DEBUG) << "(Hillclimbing optimizer)";
      transitmapper::optim::HillClimbOptimizer hillcOptim(&cfg, &scorer);
      hillcOptim.optimize(&g);
    } else if (cfg.optimMethod == "anneal") {
      LOG(DEBUG) << "(Simulated annealing optimizer)";
      transitmapper::optim::SimulatedAnnealingOptimizer annealOptim(&cfg,
                                                                    &scorer);
      annealOptim.optimize(&g);
    }

    LOG(INFO) << "(stats) Total graph score AFTER optim is -- "
              << scorer.getScore() << " -- (excl. unavoidable crossings!)";
    LOG(INFO) << "(stats)   Per node graph score: "
              << scorer.getScore() / g.getNodes()->size();
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
    transitmapper::output::SvgOutput svgOut(&o, &cfg, &scorer);
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
      transitmapper::output::SvgOutput svgOut(&o, &cfg, &scorer);
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
      transitmapper::output::SvgOutput svgOut(&o, &cfg, &scorer);
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
      transitmapper::output::SvgOutput svgOut(&o, &cfg, &scorer);
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
