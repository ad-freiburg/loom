// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "./builder/Builder.h"
#include "./config/ConfigReader.h"
#include "./config/GraphBuilderConfig.h"
#include "./graph/Graph.h"
#include "./output/JsonOutput.h"
#include "gtfsparser/Parser.h"
#include "gtfsparser/gtfs/Service.h"
#include "pbutil/log/Log.h"

using namespace graphbuilder;
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
    parser.parse(&feed, cfg.inputFeedPath);

    graphbuilder::graph::Graph g("shinygraph", cfg.projectionString);
    Builder b(&cfg);

    b.consume(feed, &g);
    b.simplify(&g);

    // first pass, with strict distance values (clearing things up first)
    b.createTopologicalNodes(&g, false);
    b.removeEdgeArtifacts(&g);

    while (b.createTopologicalNodes(&g, true)) {
    }

    b.removeEdgeArtifacts(&g);
    b.removeNodeArtifacts(&g);
    b.averageNodePositions(&g);

    JsonOutput out(&cfg);
    out.print(g);
  }

  return (0);
}
