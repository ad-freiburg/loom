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

  if (cfg.randomColors) tg.fillMissingColors();

  // snap orphan stations
  tg.snapOrphanStations();

  size_t numNdsBef = tg.getNds().size();
  size_t numEdgsBef = 0;

  double lenBef = 0, lenAfter = 0;

  if (cfg.outputStats) {
    for (const auto& nd : tg.getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        numEdgsBef++;
        lenBef += e->pl().getPolyline().getLength();
      }
    }
  }

  size_t statFr = mc.freeze();

  si.init();

  mc.averageNodePositions();


  // does preserve existing turn restrictions
  mc.removeNodeArtifacts(false);

  mc.cleanUpGeoms();

  // only remove the artifacts after the restriction inferrer has been
  // initialized, as these operations do not guarantee that the restrictions
  // are preserved!


  ri.init();
  size_t restrFr = mc.freeze();

  // util::geo::output::GeoGraphJsonOutput gout;
  // gout.printLatLng(tg, std::cout);
  // exit(0);

  mc.removeEdgeArtifacts();

  T_START(construction);
  size_t iters = 0;
  iters += mc.collapseShrdSegs(10);
  iters += mc.collapseShrdSegs(cfg.maxAggrDistance);
  double constrT = T_STOP(construction);

  mc.removeNodeArtifacts(false);

  double avgMergedEdgs = 0;
  size_t maxMergedEdgs = 0;
  if (cfg.outputStats) {
    size_t c = 0;
    const auto& origEdgs = mc.freezeTrack(restrFr);
    for (const auto& nd : tg.getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        size_t cur = origEdgs.at(e).size();
        if (cur > maxMergedEdgs) maxMergedEdgs = cur;
        avgMergedEdgs += cur;
        c++;
      }
    }
    avgMergedEdgs /= c;
  }

  mc.reconstructIntersections();


  // infer restrictions
  T_START(restrInf);
  if (!cfg.noInferRestrs) ri.infer(mc.freezeTrack(restrFr));
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

  // remove orphan lines again
  mc.removeOrphanLines();

  size_t numEdgsAfter = 0;

  if (cfg.outputStats) {
    for (const auto& nd : tg.getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        lenAfter += e->pl().getPolyline().getLength();
        numEdgsAfter++;
      }
    }
  }

  int numComps = -1;

  if (cfg.writeComponents || !cfg.componentsPath.empty()) {
    util::geo::output::GeoGraphJsonOutput out;
    const auto& graphs = tg.distConnectedComponents(10000, cfg.writeComponents);

    numComps = graphs.size();

    for (size_t comp = 0; comp < graphs.size(); comp++) {

      std::ofstream f;
      f.open(cfg.componentsPath + "/component-" + std::to_string(comp) + ".json");

      out.printLatLng(graphs[comp], f);
    }
  }


  // output
  util::geo::output::GeoGraphJsonOutput out;
  if (cfg.outputStats) {
    util::json::Dict jsonStats = {
        {"statistics",
         util::json::Dict{
             {"num_edgs_in", numEdgsBef},
             {"num_nds_in", numNdsBef},
             {"num_edgs_out", numEdgsAfter},
             {"num_nds_out", tg.getNds().size()},
             {"num_components", numComps},
             {"time_const", constrT},
             {"iters", iters},
             {"time_const", constrT},
             {"time_restr_inf", restrT},
             {"time_station_insert", stationT},
             {"len_before", lenBef},
             {"num_restrs", tg.numConnExcs()},
             {"avg_merged_edgs", avgMergedEdgs},
             {"max_merged_edgs", maxMergedEdgs},
             {"len_after", lenAfter},
         }}};
    out.printLatLng(tg, std::cout, jsonStats);
  } else {
    out.printLatLng(tg, std::cout);
  }

  return (0);
}
