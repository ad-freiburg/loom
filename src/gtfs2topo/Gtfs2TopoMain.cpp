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
#include "util/geo/output/GeoGraphJsonOutput.h"
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

    std::cerr << "Consuming... " << std::endl;
    b.consume(feed, &g);
    std::cerr << "Simplifying... " << std::endl;
    b.simplify(&g);

    // first pass, with strict distance values (clearing things up first)
    std::cerr << "Topo 1... " << std::endl;
    b.createTopologicalNodes(&g, false);
    std::cerr << "Removing artifacts... " << std::endl;
    b.removeEdgeArtifacts(&g);

    std::cerr << "Topo 2... " << std::endl;
    while (b.createTopologicalNodes(&g, true)) {
    }

    std::cerr << "Removing artifacts... " << std::endl;
    b.removeEdgeArtifacts(&g);
    b.removeNodeArtifacts(&g);
    b.averageNodePositions(&g);

    util::geo::output::GeoGraphJsonOutput out;
    out.print(g, std::cout);
  }

  return (0);
}
