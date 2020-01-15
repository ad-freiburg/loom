// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "shared/linegraph/LineGraph.h"
#include "topo/mapconstructor/MapConstructor.h"
#include "topo/statinserter/StatInserter.h"
#include "topo/config/ConfigReader.h"
#include "topo/config/TopoConfig.h"
#include "topo/restr/RestrInferrer.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  topo::config::TopoConfig cfg;
  shared::linegraph::LineGraph tg;
  topo::restr::RestrInferrer ri(&cfg, &tg);
  topo::MapConstructor mc(&cfg, &tg);
  topo::StatInserter si(&cfg, &tg);

  // read config
  topo::config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  // read input graph
  tg.readFromJson(&(std::cin), -1);

  size_t statFr = mc.freeze();
  si.init();

  mc.averageNodePositions();
  mc.cleanUpGeoms();

  mc.removeNodeArtifacts();
  mc.removeEdgeArtifacts();

  // init restriction inferrer
  ri.init();
  size_t restrFr = mc.freeze();

  // first run, with 0 perc of line width, and offset of 5
  mc.collapseShrdSegs(5.0);

  double step = cfg.maxAggrDistance;

  for (double d = cfg.maxAggrDistance; d <= (cfg.maxAggrDistance * 15);
       d += step) {
    std::cerr << d << std::endl;
    while (mc.collapseShrdSegs(d)) {
      mc.removeNodeArtifacts();
      mc.removeEdgeArtifacts();
    };
  }

  mc.removeNodeArtifacts();
  mc.averageNodePositions();

  // infer restrictions
  ri.infer(mc.freezeTrack(restrFr));

  // insert stations
  si.insertStations(mc.freezeTrack(statFr));

  // output
  util::geo::output::GeoGraphJsonOutput out;
  out.print(tg, std::cout);

  return (0);
}
