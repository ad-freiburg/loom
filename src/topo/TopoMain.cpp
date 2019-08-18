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

  // TODO: read graph

  // Builder b(&cfg);

  // b.simplify(&g);

  // // first pass, with strict distance values (clearing things up first)
  // b.createTopologicalNodes(&g, false);
  // b.removeEdgeArtifacts(&g);

  // while (b.createTopologicalNodes(&g, true)) {
  // }

  // b.removeEdgeArtifacts(&g);
  // b.removeNodeArtifacts(&g);
  // b.averageNodePositions(&g);

  util::geo::output::GeoGraphJsonOutput out;
  out.print(tg, std::cout);

  return (0);
}
