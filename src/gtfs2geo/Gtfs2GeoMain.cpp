// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "gtfs2geo/builder/Builder.h"
#include "gtfs2geo/config/ConfigReader.h"
#include "gtfs2geo/config/GraphBuilderConfig.h"
#include "gtfs2geo/graph/BuildGraph.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/gtfs/Service.h"
#include "util/log/Log.h"
#include "gtfs2geo/graph/EdgePL.h"
#include "gtfs2geo/graph/NodePL.h"

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

  // parse an example feed
  ad::cppgtfs::Parser parser;
  ad::cppgtfs::gtfs::Feed feed;

  if (!cfg.inputFeedPath.empty()) {
    parser.parse(&feed, cfg.inputFeedPath);
    gtfs2geo::graph::BuildGraph g;
    Builder b(&cfg);

    b.consume(feed, &g);
    b.simplify(&g);

    util::geo::output::GeoGraphJsonOutput out;
    out.print(g, std::cout);
  }

  return (0);
}
