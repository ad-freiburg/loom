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

  double lenBef = 0, lenAfter = 0;

  if (cfg.outputStats) {
    for (const auto& nd : *tg.getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        lenBef += e->pl().getPolyline().getLength();
      }
    }
  }

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

  T_START(construction);
  size_t iters = 0;
  iters += mc.collapseShrdSegs(10);
  iters += mc.collapseShrdSegs(cfg.maxAggrDistance);
  double constrT = T_STOP(construction);

  mc.removeNodeArtifacts(false);

  // infer restrictions
  T_START(restrInf);
  ri.infer(mc.freezeTrack(restrFr));
  double restrT = T_STOP(restrInf);

  // insert stations
  T_START(stationIns);
  si.insertStations(mc.freezeTrack(statFr));
  double stationT = T_STOP(stationIns);

  // remove orphan lines, which may be introduced by another station
  // placement
  mc.removeOrphanLines();

  mc.removeNodeArtifacts(true);

  mc.reconstructIntersections();

  if (cfg.outputStats) {
    for (const auto& nd : *tg.getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        lenAfter += e->pl().getPolyline().getLength();
      }
    }
  }

  // output
  util::geo::output::GeoGraphJsonOutput out;
  if (cfg.outputStats) {
    util::json::Dict jsonStats = {
        {"statistics",
         util::json::Dict{
             {"iters", iters},
             {"time_const", constrT},
             {"time_restr_inf", restrT},
             {"time_station_insert", stationT},
             {"len_before", lenBef},
             {"num_restrs", tg.numConnExcs()},
             {"len_after", lenAfter},
         }}};
    out.print(tg, std::cout, jsonStats);
  } else {
    out.print(tg, std::cout);
  }

  return (0);
}
