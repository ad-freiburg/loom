// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/builder/Builder.h"
#include "topo/config/ConfigReader.h"
#include "topo/config/TopoConfig.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  topo::config::TopoConfig cfg;

  topo::config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  shared::transitgraph::TransitGraph tg;
  tg.readFromJson(&(std::cin));

  std::cerr << "Topoying graph with " << tg.getNds()->size() << " nodes..."
            << std::endl;

  topo::Builder b(&cfg);


  b.initEdges(&tg);
  b.collectStations(&tg);

  b.averageNodePositions(&tg);
  b.cleanUpGeoms(&tg);

  b.removeNodeArtifacts(&tg);
  b.removeEdgeArtifacts(&tg);

  // first pass, with strict distance values (clearing things up first)
  std::cerr << "Creating topological nodes (first round)..." << std::endl;

  // first run, with 0 perc of line width, and offset of 5
  b.createTopologicalNodes(&tg, 5.0);

  double step = cfg.maxAggrDistance;

  for (double d = cfg.maxAggrDistance; d <= (cfg.maxAggrDistance * 15); d += step) {
    std::cerr << d  << std::endl;
    while (b.createTopologicalNodes(&tg, d)) {
      b.removeNodeArtifacts(&tg);
      b.removeEdgeArtifacts(&tg);
    };
  }

  std::cerr << tg.getNds()->size() << " nodes..." << std::endl;
  std::cerr << "Removing node artifacts..." << std::endl;
  b.removeNodeArtifacts(&tg);
  std::cerr << tg.getNds()->size() << " nodes..." << std::endl;
  std::cerr << "Averaging node positions..." << std::endl;

  b.averageNodePositions(&tg);
  b.insertStations(&tg);

  util::geo::output::GeoGraphJsonOutput out;
  out.print(tg, std::cout);

  return (0);
}
