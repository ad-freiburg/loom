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
#include "topo/config/ConfigReader.h"
#include "topo/config/TopoConfig.h"
#include "topo/mapconstructor/MapConstructor.h"
#include "topo/restr/RestrInferrer.h"
#include "topo/statinserter/StatInserter.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  topo::config::TopoConfig cfg;

  size_t iters = 0;
  double constrT = 0;
  double restrT = 0;
  double stationT = 0;

  shared::linegraph::LineGraph lg;
  // read config
  topo::config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  // read input graph
  lg.readFromJson(&(std::cin));

  if (cfg.randomColors) lg.fillMissingColors();

  // snap orphan stations
  lg.snapOrphanStations();

  size_t numNdsBef = 0;
  size_t numEdgsBef = 0;
  double lenBef = 0, lenAfter = 0;
  size_t totMergedEdgs = 0;
  size_t totSupportGraphEdgs = 0;
  size_t maxMergedEdgs = 0;

  if (cfg.outputStats) {
    if (cfg.aggregateStats) {
      const auto& props = lg.getGraphProps();
      if (props.count("statistics")) {
        const auto& stats =
            props.at("statistics").get<nlohmann::json::object_t>();
        if (stats.count("num_nds_in"))
          numNdsBef = stats.at("num_nds_in").get<size_t>();
        if (stats.count("num_edgs_in"))
          numEdgsBef = stats.at("num_edgs_in").get<size_t>();
        if (stats.count("len_before"))
          lenBef = stats.at("len_before").get<double>();
        if (stats.count("iters")) iters = stats.at("iters").get<size_t>();
        if (stats.count("time_const"))
          constrT = stats.at("time_const").get<double>();
        if (stats.count("time_restr_inf"))
          restrT = stats.at("time_restr_inf").get<double>();
        if (stats.count("time_station_insert"))
          stationT = stats.at("time_station_insert").get<double>();
        if (stats.count("max_merged_edgs"))
          maxMergedEdgs = stats.at("max_merged_edgs").get<size_t>();
        if (stats.count("tot_merged_edgs"))
          totMergedEdgs = stats.at("tot_merged_edgs").get<size_t>();
        if (stats.count("tot_support_graph_edgs"))
          totSupportGraphEdgs =
              stats.at("tot_support_graph_edgs").get<size_t>();
      }
    } else {
      numNdsBef = lg.getNds().size();
      for (const auto& nd : lg.getNds()) {
        for (const auto& e : nd->getAdjList()) {
          if (e->getFrom() != nd) continue;
          numEdgsBef++;
          lenBef += e->pl().getPolyline().getLength();
        }
      }
    }
  }

  size_t numEdgsAfter = 0;
  size_t numNdsAfter = 0;
  size_t numStationsAfter = 0;

  size_t numConExc = 0;

  lg.removeDeg1Nodes();

  LOGTO(DEBUG, std::cerr) << "Computing components...";
  auto graphs = lg.distConnectedComponents(cfg.connectedCompDist, false);

  LOGTO(DEBUG, std::cerr) << "Broke up input into " << graphs.size()
                          << " components (including single-node components)";

  std::vector<LineGraph*> resultGraphs;

  size_t compI = 0;

  // TODO: parallelize this (would increase memory usage significantly)?
  for (auto& tg : graphs) {
    LOGTO(DEBUG, std::cerr) << "@ Component" << compI++ << " components";

    topo::restr::RestrInferrer ri(&cfg, &tg);
    topo::MapConstructor mc(&cfg, &tg);
    topo::StatInserter si(&cfg, &tg);

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

    mc.removeEdgeArtifacts();

    T_START(construction);
    iters += mc.collapseShrdSegs(10, 50, cfg.segmentLength);
    iters += mc.collapseShrdSegs(cfg.maxAggrDistance, 50, cfg.segmentLength);
    constrT += T_STOP(construction);

    mc.removeNodeArtifacts(false);

    if (cfg.outputStats) {
      const auto& origEdgs = mc.freezeTrack(restrFr);
      for (const auto& nd : tg.getNds()) {
        for (const auto& e : nd->getAdjList()) {
          if (e->getFrom() != nd) continue;
          size_t cur = origEdgs.at(e).size();
          if (cur > maxMergedEdgs) maxMergedEdgs = cur;
          totMergedEdgs += cur;
          totSupportGraphEdgs++;
        }
      }
    }

    mc.reconstructIntersections();

    // infer restrictions
    T_START(restrInf);
    if (!cfg.noInferRestrs) ri.infer(mc.freezeTrack(restrFr));
    restrT += T_STOP(restrInf);

    // insert stations
    T_START(stationIns);
    si.insertStations(mc.freezeTrack(statFr));
    stationT += T_STOP(stationIns);

    // remove orphan lines, which may be introduced by another station
    // placement
    mc.removeOrphanLines();

    mc.removeNodeArtifacts(true);

    mc.reconstructIntersections();

    // remove orphan lines again
    mc.removeOrphanLines();

    if (cfg.outputStats) {
      for (const auto& nd : tg.getNds()) {
        numNdsAfter++;
        if (nd->pl().stops().size()) numStationsAfter++;
        for (const auto& e : nd->getAdjList()) {
          if (e->getFrom() != nd) continue;
          lenAfter += e->pl().getPolyline().getLength();
          numEdgsAfter++;
        }
      }
    }

    numConExc += tg.numConnExcs();

    if (cfg.smooth > 0) tg.smooth(cfg.smooth);

    resultGraphs.push_back(&tg);
  }

  int numComps = 0;

  size_t offset = 0;

  for (auto& tg : resultGraphs) {
    if (tg->getNds().size() == 0) continue;
    if (cfg.writeComponents || !cfg.componentsPath.empty()) {
      util::geo::output::GeoGraphJsonOutput out;

      size_t locOffset = offset;
      const auto& graphs = tg->distConnectedComponents(
          cfg.connectedCompDist, cfg.writeComponents, &offset);

      numComps += graphs.size();

      for (size_t comp = 0; comp < graphs.size(); comp++) {
        std::ofstream f;
        f.open(cfg.componentsPath + "/component-" +
               std::to_string(locOffset + comp) + ".json");

        out.printLatLng(graphs[comp], f);
      }
    }
  }

  // output
  util::geo::output::GeoGraphJsonOutput gout;
  if (cfg.outputStats) {
    util::json::Dict jsonStats = {
        {"statistics",
         util::json::Dict{
             {"num_edgs_in", numEdgsBef},
             {"num_nds_in", numNdsBef},
             {"num_edgs_out", numEdgsAfter},
             {"num_nds_out", numNdsAfter},
             {"num_stations_out", numStationsAfter},
             {"num_components", numComps},
             {"time_const", constrT},
             {"iters", iters},
             {"time_const", constrT},
             {"time_restr_inf", restrT},
             {"time_station_insert", stationT},
             {"len_before", lenBef},
             {"num_restrs", numConExc},
             {"avg_merged_edgs", (static_cast<double>(totMergedEdgs) /
                                  static_cast<double>(totSupportGraphEdgs))},
             {"max_merged_edgs", maxMergedEdgs},
             {"len_after", lenAfter},
             {"tot_merged_edgs", totMergedEdgs},
             {"tot_support_graph_edgs", totSupportGraphEdgs},
         }}};

    util::geo::output::GeoJsonOutput out(std::cout, jsonStats);
    for (auto gg : resultGraphs) {
      gout.printLatLng(*gg, &out);
    }
    out.flush();
  } else {
    util::geo::output::GeoJsonOutput out(std::cout);
    for (auto gg : resultGraphs) {
      gout.printLatLng(*gg, &out);
    }
    out.flush();
  }

  return (0);
}
