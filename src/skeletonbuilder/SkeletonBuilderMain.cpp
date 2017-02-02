// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <stdio.h>
#include "log/Log.h"
#include "gtfsparser/Parser.h"
#include "transitmap/graph/TransitGraph.h"
#include "./builder/Builder.h"
#include "transitmap/output/OgrOutput.h"
#include "gtfsparser/gtfs/Service.h"
#include "transitmap/config/ConfigReader.cpp"
#include "transitmap/config/TransitMapConfig.h"

using namespace skeletonbuilder;
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

  // parse an example feed
  gtfsparser::Parser parser;
  gtfsparser::gtfs::Feed feed;

  if (!cfg.inputFeedPath.empty()) {
    LOG(INFO) << "reading feed at " << cfg.inputFeedPath << std::endl;
    parser.parse(&feed, cfg.inputFeedPath);

    graph::TransitGraph g("shinygraph", cfg.projectionString);
    Builder b(&cfg);

    LOG(INFO) << "Building graph..." << std::endl;
    b.consume(feed, &g);

    LOG(INFO) << "Simplyfing..." << std::endl;
    b.simplify(&g);

    LOG(INFO) << "Building topological nodes..." << std::endl;
    while (b.createTopologicalNodes(&g)) {}

    LOG(INFO) << "Averaging node positions" << std::endl;
    b.averageNodePositions(&g);

    b.removeArtifacts(&g);

    std::string path = cfg.outputPath;
    LOG(INFO) << "Outputting to OGR " << path << " ..."
      << std::endl;
    output::OgrOutput ogrOut(path, &cfg);
    ogrOut.print(g);

  }

  return(0);
}
