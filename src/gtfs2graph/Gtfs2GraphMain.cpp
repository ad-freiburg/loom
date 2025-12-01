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
#include "gtfs2graph/builder/Builder.h"
#include "gtfs2graph/config/ConfigReader.h"
#include "gtfs2graph/config/GraphBuilderConfig.h"
#include "gtfs2graph/graph/BuildGraph.h"
#include "gtfs2graph/graph/EdgePL.h"
#include "gtfs2graph/graph/NodePL.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using namespace gtfs2graph;
using util::ERROR;
using util::DEBUG;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  // parse an example feed
  ad::cppgtfs::gtfs::Feed feed;

  if (!cfg.inputFeedPath.empty()) {
    try {
      ad::cppgtfs::Parser parser(cfg.inputFeedPath);
      parser.parse(&feed);
    } catch (const ad::cppgtfs::ParserException& ex) {
      LOG(ERROR) << "Could not parse input GTFS feed, reason was:";
      std::cerr << ex.what() << std::endl;
      exit(1);
    }
    gtfs2graph::graph::BuildGraph g;
    Builder b(&cfg);

    b.consume(feed, &g);

    b.simplify(&g);

    util::geo::output::GeoGraphJsonOutput out;
    out.printLatLng(g, std::cout);
  }

  return 0;
}
