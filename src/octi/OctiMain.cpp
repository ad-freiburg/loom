// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi <brosi@cs.uni-freiburg.de>

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include "json/json.hpp"
#include "octi/Octilinearizer.h"
#include "octi/Enlarger.h"
#include "octi/basegraph/BaseGraph.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/config/ConfigReader.h"
#include "shared/linegraph/LineGraph.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/BiDijkstra.h"
#include "util/json/Writer.h"
#include "util/log/Log.h"
#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_num_procs() 1
#endif

using std::string;
using namespace octi;

using octi::Octilinearizer;
using octi::Enlarger;
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
const CombNode* getCenterNd(const CombGraph* cg) {
  const CombNode* ret = 0;
  for (auto nd : cg->getNds()) {
    if (!ret ||
        LineGraph::getLDeg(nd->pl().getParent()) >
            LineGraph::getLDeg(ret->pl().getParent())) {
      ret = nd;
    }
  }

  return ret;
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
    LOGTO(DEBUG, std::cerr) << "Reading obstacle file... ";
    cfg.obstacles = readObstacleFile(cfg.obstaclePath);
    LOGTO(DEBUG, std::cerr)
        << "Done. (" << cfg.obstacles.size() << " obstacles)";
  }

  LOGTO(DEBUG, std::cerr) << "Reading graph file... ";
  T_START(read);
  LineGraph tg;
  BaseGraph* gg;
  if (cfg.fromDot)
    tg.readFromDot(&(std::cin), 0);
  else
    tg.readFromJson(&(std::cin), 0);
  LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(read) << "ms)";

  LOGTO(DEBUG, std::cerr) << "Planarizing graph... ";
  T_START(planarize);
  tg.topologizeIsects();
  LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(planarize) << "ms)";

  double avgDist = avgStatDist(tg);
  LOGTO(DEBUG, std::cerr) << "Average adj. node distance is " << avgDist;

  // BaseGraphType graphType = BaseGraphType::ORTHORADIAL;
  // BaseGraphType graphType = BaseGraphType::OCTIHANANGRID;
  // BaseGraphType graphType = BaseGraphType::OCTIQUADTREE;
  BaseGraphType graphType = BaseGraphType::OCTIGRID;
  // BaseGraphType graphType = BaseGraphType::GRID;

  Octilinearizer oct(graphType);
  LineGraph res;

  double gridSize;

  if (util::trim(cfg.gridSize).back() == '%') {
    double perc = atof(cfg.gridSize.c_str()) / 100;
    gridSize = avgDist * perc;
    LOGTO(DEBUG, std::cerr)
        << "Grid size " << gridSize << " (" << perc * 100 << "%)";
  } else {
    gridSize = atof(cfg.gridSize.c_str());
    LOGTO(DEBUG, std::cerr) << "Grid size " << gridSize;
  }

  tg.contractEdges(gridSize / 2);

  auto box = tg.getBBox();

  // TODO!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  // change according to graphtype
  tg.splitNodes(4);
  CombGraph cg(&tg, cfg.deg2Heur);

  // local enlargement
  LOGTO(DEBUG, std::cerr) << "Locally enlarging graph... ";
  Enlarger e;
  e.enlarge(cg, 100);
  avgDist = avgStatDist(tg);
  LOGTO(DEBUG, std::cerr)
      << "Average adj. node distance after local enlargement is " << avgDist;

  box = util::geo::pad(box, gridSize + 1);

  if (graphType == ORTHORADIAL) {
    auto centerNd = getCenterNd(&cg);

    std::cerr << "Center node is "
              << centerNd->pl().getParent()->pl().toString() << std::endl;

    auto cgCtr = *centerNd->pl().getGeom();
    auto newBox = util::geo::DBox();

    newBox = extendBox(box, newBox);
    newBox = extendBox(rotate(convexHull(box), 180, cgCtr), newBox);
    box = newBox;
  }

  Score sc;
  double time = 0;

  if (cfg.optMode == "ilp") {
    T_START(octi);
    sc = oct.drawILP(cg, box, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                     cfg.maxGrDist, cfg.ilpNoSolve, cfg.enfGeoPen,
                     cfg.ilpTimeLimit, cfg.ilpSolver, cfg.ilpPath);
    time = T_STOP(octi);
    LOGTO(DEBUG, std::cerr)
        << "Octilinearized using ILP in " << time << " ms, score " << sc.full;
  } else if ((cfg.optMode == "heur")) {
    T_START(octi);
    try {
      sc = oct.draw(cg, box, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                    cfg.maxGrDist, cfg.restrLocSearch,
                    cfg.enfGeoPen, cfg.obstacles, cfg.abortAfter);
      time = T_STOP(octi);
    } catch (const NoEmbeddingFoundExc& exc) {
      LOG(ERROR) << exc.what();
      exit(1);
    }
    LOGTO(DEBUG, std::cerr) << "Octilinearized using heur approach in " << time
                            << " ms, score " << sc.full;
  }

  size_t maxRss = util::getPeakRSS();

  // translate score to JSON
  util::json::Dict jsonScore{
      {"scores", util::json::Dict{{"total_score", sc.full},
                                  {"density-score", sc.dense},
                                  {"bend-score", sc.bend},
                                  {"hop-score", sc.hop},
                                  {"move-score", sc.move}}},
      {"pens",
       util::json::Dict{
           {"density-pen", cfg.pens.densityPen},
           {"diag-pen", cfg.pens.diagonalPen},
           {"hori-pen", cfg.pens.horizontalPen},
           {"vert-pen", cfg.pens.verticalPen},
           {"180-turn-pen", cfg.pens.p_0},
           {"135-turn-pen", cfg.pens.p_135},
           {"90-turn-pen", cfg.pens.p_90},
           {"45-turn-pen", cfg.pens.p_45},
       }},
      {"misc", util::json::Dict{{"method", cfg.optMode},
                                {"deg2heur", cfg.deg2Heur},
                                {"max-grid-dist", cfg.maxGrDist}}},
      {"time_ms", time},
      {"procs", omp_get_num_procs()},
      {"peak_memory", util::readableSize(maxRss)},
      {"peak_memory_bytes", util::json::Int(maxRss)},
      {"timestamp", util::json::Int(std::time(0))}};

  if (cfg.printMode == "gridgraph") {
    out.print(*gg, std::cout, jsonScore);
  } else {
    out.print(res, std::cout, jsonScore);
  }

  return (0);
}
