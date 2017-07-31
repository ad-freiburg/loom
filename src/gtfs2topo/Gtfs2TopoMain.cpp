// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "gtfs2topo/builder/Builder.h"
#include "gtfs2topo/config/ConfigReader.h"
#include "gtfs2topo/config/GraphBuilderConfig.h"
#include "gtfs2topo/graph/BuildGraph.h"
#include "gtfs2topo/output/GeoJsonOutput.h"
#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/gtfs/Service.h"
#include "util/log/Log.h"
#include "gtfs2topo/graph/EdgePL.h"
#include "gtfs2topo/graph/NodePL.h"

using namespace gtfs2topo;
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
  ad::cppgtfs::Parser parser;
  ad::cppgtfs::gtfs::Feed feed;

  if (!cfg.inputFeedPath.empty()) {
    parser.parse(&feed, cfg.inputFeedPath);

    gtfs2topo::graph::BuildGraph g;
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

    GeoJsonOutput out(&cfg);
    out.print(g);
  }

  return (0);
}
