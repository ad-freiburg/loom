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

using namespace topo;
using std::string;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  config::TopoConfig cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  shared::transitgraph::TransitGraph tg;
  tg.readFromJson(&(std::cin));

  Builder b(&cfg);

  // first pass, with strict distance values (clearing things up first)
  b.createTopologicalNodes(&tg, false);
  b.removeEdgeArtifacts(&tg);

  while (b.createTopologicalNodes(&tg, true));

  b.removeEdgeArtifacts(&tg);
  b.removeNodeArtifacts(&tg);
  b.averageNodePositions(&tg);

  util::geo::output::GeoGraphJsonOutput out;
  out.print(tg, std::cout);

  return (0);
}
