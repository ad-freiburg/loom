// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include "json/json.hpp"
#include "octi/Octilinearizer.h"
#include "octi/basegraph/BaseGraph.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/config/ConfigReader.h"
#include "shared/linegraph/LineGraph.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/BiDijkstra.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

using octi::Octilinearizer;
using octi::basegraph::BaseGraph;
using util::geo::dist;
using util::geo::DPolygon;

// _____________________________________________________________________________
double avgStatDist(const LineGraph& g) {
  double avg = 0;
  size_t i = 0;
  for (const auto nd : g.getNds()) {
    if (nd->getDeg() == 0) continue;
    i++;
    double loc = 0;
    for (const auto edg : nd->getAdjList()) {
      loc += dist(*nd->pl().getGeom(), *edg->getOtherNd(nd)->pl().getGeom());
    }
    avg += loc / nd->getAdjList().size();
  }
  avg /= i++;
  return avg;
}

// _____________________________________________________________________________
std::vector<DPolygon> readObstacleFile(const std::string& p) {
  std::vector<DPolygon> ret;
  std::ifstream s;
  s.open(p);
  nlohmann::json j;
  s >> j;

  if (j["type"] == "FeatureCollection") {
    for (auto feature : j["features"]) {
      auto geom = feature["geometry"];
      if (geom["type"] == "Polygon") {
        std::vector<std::vector<double>> coords = geom["coordinates"][0];
        util::geo::Line<double> l;
        for (auto coord : coords) {
          l.push_back({coord[0], coord[1]});
        }
        ret.push_back(l);
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  util::geo::output::GeoGraphJsonOutput out;

  if (cfg.obstaclePath.size()) {
    LOG(INFO, std::cerr) << "Reading obstacle file... ";
    cfg.obstacles = readObstacleFile(cfg.obstaclePath);
    LOG(INFO, std::cerr) << "Done. (" << cfg.obstacles.size() << " obstacles)";
  }

  LOG(INFO, std::cerr) << "Reading graph file... ";
  T_START(read);
  LineGraph tg;
  BaseGraph* gg;
  if (cfg.fromDot)
    tg.readFromDot(&(std::cin), 0);
  else
    tg.readFromJson(&(std::cin), 0);
  LOG(INFO, std::cerr) << "Done. (" << T_STOP(read) << "ms)";

  LOG(INFO, std::cerr) << "Planarize graph... ";
  T_START(planarize);
  tg.topologizeIsects();
  LOG(INFO, std::cerr) << "Done. (" << T_STOP(planarize) << "ms)";

  double avgDist = avgStatDist(tg);
  LOG(INFO, std::cerr) << "Average adj. node distance is " << avgDist;

  // BaseGraphType graphType = BaseGraphType::ORTHORADIAL;
  BaseGraphType graphType = BaseGraphType::OCTIHANANGRID;
  // BaseGraphType graphType = BaseGraphType::OCTIQUADTREE;
  // BaseGraphType graphType = BaseGraphType::OCTIGRID;
  // BaseGraphType graphType = BaseGraphType::GRID;

  Octilinearizer oct(graphType);
  LineGraph res;

  double gridSize;

  if (util::trim(cfg.gridSize).back() == '%') {
    double perc = atof(cfg.gridSize.c_str()) / 100;
    gridSize = avgDist * perc;
    LOG(INFO, std::cerr) << "Grid size " << gridSize << " (" << perc * 100
                         << "%)";
  } else {
    gridSize = atof(cfg.gridSize.c_str());
    LOG(INFO, std::cerr) << "Grid size " << gridSize;
  }

  if (cfg.optMode == "ilp") {
    T_START(octi);
    double sc =
        oct.drawILP(&tg, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                    cfg.deg2Heur, cfg.maxGrDist, cfg.ilpNoSolve, cfg.enfGeoPen,
                    cfg.ilpTimeLimit, cfg.ilpSolver, cfg.ilpPath);
    LOG(INFO, std::cerr) << "Octilinearized using ILP in " << T_STOP(octi)
                         << " ms, score " << sc;
  } else if ((cfg.optMode == "heur")) {
    T_START(octi);
    double sc;
    try {
      sc = oct.draw(&tg, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                    cfg.deg2Heur, cfg.maxGrDist, cfg.restrLocSearch,
                    cfg.enfGeoPen, cfg.obstacles);
    } catch (const NoEmbeddingFoundExc& exc) {
      LOG(ERROR) << exc.what();
      exit(1);
    }
    LOG(INFO, std::cerr) << "Octilinearized using heur approach in "
                         << T_STOP(octi) << " ms, score " << sc;
  }

  if (cfg.printMode == "gridgraph") {
    out.print(*gg, std::cout);
  } else {
    out.print(res, std::cout);
  }

  return (0);
}
