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
  tg.readFromJson(&(std::cin), 0);

  size_t statFr = mc.freeze();
  si.init();

  mc.averageNodePositions();

  mc.cleanUpGeoms();

  // does preserve existing turn restrictions
  mc.removeNodeArtifacts(false);

  // init restriction inferrer
  ri.init();
  size_t restrFr = mc.freeze();

  // only remove the artifacts after the restriction inferrer has been
  // initialized, as these operations do not guarantee that the restrictions
  // are preserved!
  mc.removeEdgeArtifacts();

  // segment length
  mc.collapseShrdSegs(10);

  mc.collapseShrdSegs(cfg.maxAggrDistance);

  mc.removeNodeArtifacts(false);

  // infer restrictions
  ri.infer(mc.freezeTrack(restrFr));

  // insert stations
  si.insertStations(mc.freezeTrack(statFr));

  // remove orphan lines, which may be introduced by another station
  // placement
  mc.removeOrphanLines();

  mc.removeNodeArtifacts(true);

  mc.reconstructIntersections();

  // output
  util::geo::output::GeoGraphJsonOutput out;
  out.print(tg, std::cout);

  return (0);
}
