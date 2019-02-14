// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/gtfs/Service.h"
#include "gtfs2geo/builder/Builder.h"
#include "gtfs2geo/config/ConfigReader.h"
#include "gtfs2geo/config/GraphBuilderConfig.h"
#include "gtfs2geo/graph/BuildGraph.h"
#include "gtfs2geo/graph/EdgePL.h"
#include "gtfs2geo/graph/NodePL.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using namespace gtfs2geo;
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

  if (!cfg.inputFeedPath.empty()) {
    gtfs2geo::graph::BuildGraph g;
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

    util::geo::output::GeoGraphJsonOutput out;
    out.print(g, std::cout);
  }

  return (0);
}
