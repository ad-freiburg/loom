// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi <brosi@cs.uni-freiburg.de>

#include <stdio.h>
#include <unistd.h>

#include <fstream>
#include <iostream>
#include <set>

#include "3rdparty/json.hpp"
#include "octi/Enlarger.h"
#include "octi/Octilinearizer.h"
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

using octi::Enlarger;
using octi::Octilinearizer;
using octi::basegraph::BaseGraph;
using util::geo::dist;
using util::geo::DPolygon;

struct TotalScore {
  Score score;
  octi::ilp::ILPStats ilpstats;

  size_t gridgraphNumNds = 0;
  size_t gridgraphNumEdgs = 0;
  size_t combgraphNumNds = 0;
  size_t combgraphNumEdgs = 0;
  size_t inputgraphNumNds = 0;
  size_t inputgraphNumEdgs = 0;
  size_t inputgraphMaxDeg = 0;
  size_t numNoEmbeddingFound = 0;
  double timeMs = 0;
};

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
    if (!ret || LineGraph::getLDeg(nd->pl().getParent()) >
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
void drawComp(LineGraph& tg, double avgDist, util::json::Array& jsonScores,
              std::vector<LineGraph*>& resultGraphs,
              std::vector<BaseGraph*>& resultGridGraphs, TotalScore& totScore,
              const config::Config& cfg) {
  Drawing d;

  Octilinearizer oct(cfg.baseGraphType);
  LineGraph* res = new LineGraph();
  BaseGraph* gg;

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

  // contract degree 2 nodes without any significance (no station, no
  // exception, no change in lines
  tg.contractStrayNds();

  // heuristic: contract all edges shorter than half the grid size
  tg.contractEdges(gridSize / 2);

  auto box = tg.getBBox();

  // split nodes that have a larger degree than the max degree of the grid
  // graph to allow drawing
  tg.splitNodes(oct.maxNodeDeg());

  CombGraph cg(&tg, cfg.deg2Heur);
  box = util::geo::pad(box, gridSize + 1);

  if (cfg.baseGraphType == octi::basegraph::BaseGraphType::ORTHORADIAL ||
      cfg.baseGraphType == octi::basegraph::BaseGraphType::PSEUDOORTHORADIAL) {
    auto centerNd = getCenterNd(&cg);

    LOGTO(DEBUG, std::cerr) << "Orthoradial center node is "
                            << centerNd->pl().getParent()->pl().toString();

    auto cgCtr = *centerNd->pl().getGeom();
    auto newBox = util::geo::DBox();

    newBox = extendBox(box, newBox);
    newBox = extendBox(rotate(convexHull(box), 180, cgCtr), newBox);
    box = newBox;
  }

  Score sc;
  octi::ilp::ILPStats ilpstats;
  double time = 0;

  if (cfg.optMode == "ilp") {
    T_START(octi);
    sc = oct.drawILP(cg, box, res, &gg, &d, cfg.pens, gridSize, cfg.borderRad,
                     cfg.maxGrDist, cfg.orderMethod, cfg.ilpNoSolve,
                     cfg.enfGeoPen, cfg.hananIters, cfg.ilpTimeLimit,
                     cfg.ilpCacheDir, cfg.ilpCacheThreshold, cfg.ilpNumThreads,
                     &ilpstats, cfg.ilpSolver, cfg.ilpPath);
    time = T_STOP(octi);
    LOGTO(DEBUG, std::cerr)
        << "Schematized using ILP in " << time << " ms, score " << sc.full;
  } else if ((cfg.optMode == "heur")) {
    T_START(octi);
    sc = oct.draw(cg, box, res, &gg, &d, cfg.pens, gridSize, cfg.borderRad,
                  cfg.maxGrDist, cfg.orderMethod, cfg.restrLocSearch,
                  cfg.enfGeoPen, cfg.hananIters, cfg.obstacles,
                  cfg.heurLocSearchIters, cfg.abortAfter);
    time = T_STOP(octi);

    LOGTO(DEBUG, std::cerr) << "Schematized using heur approach in " << time
                            << " ms, score " << sc.full;
  }

  if (cfg.writeStats) {
    size_t maxRss = util::getPeakRSS();
    size_t numEdgs = 0;
    size_t numEdgsComb = 0;
    size_t numEdgsTg = 0;
    for (auto nd : gg->getNds()) {
      numEdgs += nd->getDeg();
    }
    for (auto nd : cg.getNds()) {
      numEdgsComb += nd->getDeg();
    }
    for (auto nd : tg.getNds()) {
      numEdgsTg += nd->getDeg();
    }

    // total score
    totScore.score = totScore.score + sc;
    totScore.ilpstats = totScore.ilpstats + ilpstats;

    totScore.gridgraphNumNds += gg->getNds().size();
    totScore.gridgraphNumEdgs += numEdgs / 2;
    totScore.combgraphNumNds += cg.getNds().size();
    totScore.combgraphNumEdgs += numEdgsComb / 2;
    totScore.inputgraphNumNds += tg.getNds().size();
    totScore.inputgraphNumEdgs += numEdgsTg / 2;
    totScore.inputgraphMaxDeg =
        std::max(totScore.inputgraphMaxDeg, tg.maxDeg());
    totScore.timeMs += time;

    // translate score to JSON
    util::json::Dict jsonScore = util::json::Dict{
        {"scores",
         util::json::Dict{{"total-score", sc.full},
                          {"topo-violations", util::json::Int(sc.violations)},
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
        {"gridgraph-size", util::json::Dict{{"nodes", gg->getNds().size()},
                                            {"edges", numEdgs / 2}}},
        {"combgraph-size", util::json::Dict{{"nodes", cg.getNds().size()},
                                            {"edges", numEdgsComb / 2}}},
        {"input-graph-size", util::json::Dict{{"nodes", tg.getNds().size()},
                                              {"edges", numEdgsTg / 2},
                                              {"max-deg", tg.maxDeg()}}},
        {"input-graph-avg-node-dist",
         avgDist * webMercDistFactor(box.getLowerLeft())},
        {"area", dist(box.getLowerRight(), box.getLowerLeft()) *
                     webMercDistFactor(box.getLowerRight()) *
                     dist(box.getLowerRight(), box.getUpperRight()) *
                     webMercDistFactor(box.getLowerRight())},
        {"misc", util::json::Dict{{"method", cfg.optMode},
                                  {"deg2heur", cfg.deg2Heur},
                                  {"max-grid-dist", cfg.maxGrDist}}},
        {"time-ms", time},
        {"iterations", sc.iters},
        {"procs", omp_get_num_procs()},
        {"peak-memory", util::readableSize(maxRss)},
        {"peak-memory-bytes", maxRss},
        {"timestamp", util::json::Int(std::time(0))}};

    if (cfg.optMode == "ilp") {
      jsonScore["ilp"] = util::json::Dict{
          {"size",
           util::json::Dict{{"rows", ilpstats.rows}, {"cols", ilpstats.cols}}},
          {"solve-time", ilpstats.time},
          {"optimal", util::json::Bool{ilpstats.optimal}}};
    }

    jsonScores.push_back(jsonScore);
  }

  resultGraphs.push_back(res);

  if (cfg.printMode == "gridgraph") {
    resultGridGraphs.push_back(gg);
  } else {
    delete gg;
  }
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
    LOGTO(DEBUG, std::cerr) << "Reading obstacle file...";
    cfg.obstacles = readObstacleFile(cfg.obstaclePath);
    LOGTO(DEBUG, std::cerr) << "Done. (" << cfg.obstacles.size() << " obst.)";
  }

  LOGTO(DEBUG, std::cerr) << "Reading graph file...";
  T_START(read);
  LineGraph lg;

  if (cfg.fromDot)
    lg.readFromDot(&(std::cin), 0);
  else
    lg.readFromJson(&(std::cin), 0);

  LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(read) << "ms)";

  LOGTO(DEBUG, std::cerr) << "Planarizing graph...";
  T_START(planarize);
  lg.topologizeIsects();
  LOGTO(DEBUG, std::cerr) << "Done. (" << T_STOP(planarize) << "ms)";

  std::vector<LineGraph> comps = lg.distConnectedComponents(10000, false);

  util::json::Array jsonScores;
  std::vector<LineGraph*> resultGraphs;
  std::vector<BaseGraph*> resultGridGraphs;

  LOGTO(DEBUG, std::cerr) << "Broke input graph into " << comps.size()
                          << " components";

  TotalScore totScore;

  for (auto& tg : comps) {
    double avgDist = avgStatDist(tg);

    double curDist = avgDist;

    size_t tries = 0;
    size_t MAX_TRIES = 10;

    LOGTO(DEBUG, std::cerr) << "Average adj. node distance is " << avgDist;

    while (tries < MAX_TRIES) {
      try {
        drawComp(tg, curDist, jsonScores, resultGraphs, resultGridGraphs,
                 totScore, cfg);

        break;
      } catch (const NoEmbeddingFoundExc& exc) {
        if (cfg.retryOnError && tries < MAX_TRIES) {
          curDist -= avgDist * 0.1;
          tries++;
          LOGTO(WARN, std::cerr) << "Retrying with grid size " << curDist;
          continue;
        }

        if (cfg.skipOnError) {
          totScore.numNoEmbeddingFound += 1;
          jsonScores.push_back(util::json::Dict());
          LOGTO(WARN, std::cerr) << exc.what();
          break;
        }

        LOG(ERROR) << exc.what();
        exit(1);
      }
    }
  }

  util::geo::output::GeoGraphJsonOutput gout;

  size_t maxRss = util::getPeakRSS();

  // translate score to JSON
  util::json::Dict totalScore = util::json::Dict{
      {"scores", util::json::Dict{{"total-score", totScore.score.full},
                                  {"topo-violations",
                                   util::json::Int(totScore.score.violations)},
                                  {"density-score", totScore.score.dense},
                                  {"bend-score", totScore.score.bend},
                                  {"hop-score", totScore.score.hop},
                                  {"move-score", totScore.score.move}}},
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

      {"gridgraph-size",
       util::json::Dict{{"nodes", totScore.gridgraphNumNds},
                        {"edges", totScore.gridgraphNumEdgs}}},
      {"combgraph-size",
       util::json::Dict{{"nodes", totScore.combgraphNumNds},
                        {"edges", totScore.combgraphNumEdgs}}},
      {"input-graph-size",
       util::json::Dict{{"nodes", totScore.inputgraphNumNds},
                        {"edges", totScore.inputgraphNumEdgs},
                        {"max-deg", totScore.inputgraphMaxDeg}}},
      {"misc", util::json::Dict{{"method", cfg.optMode},
                                {"deg2heur", cfg.deg2Heur},
                                {"max-grid-dist", cfg.maxGrDist}}},
      {"num-comps-no-embedding-found", totScore.numNoEmbeddingFound},
      {"num-comps", comps.size()},
      {"time-ms", totScore.timeMs},
      {"iterations", totScore.score.iters},
      {"procs", omp_get_num_procs()},
      {"peak-memory", util::readableSize(maxRss)},
      {"peak-memory-bytes", maxRss},
      {"timestamp", util::json::Int(std::time(0))}};

  if (cfg.optMode == "ilp") {
    totalScore["ilp"] = util::json::Dict{
        {"size", util::json::Dict{{"rows", totScore.ilpstats.rows},
                                  {"cols", totScore.ilpstats.cols}}},
        {"solve-time", totScore.ilpstats.time},
        {"optimal", util::json::Bool{totScore.ilpstats.optimal}}};
  }

  if (cfg.printMode == "gridgraph") {
    if (cfg.writeStats) {
      util::geo::output::GeoJsonOutput out(
          std::cout, util::json::Dict{{"statistics", totalScore},
                                      {"component-statistics", jsonScores}});
      for (auto gg : resultGridGraphs) {
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
  } else {
    if (cfg.writeStats) {
      util::geo::output::GeoJsonOutput out(
          std::cout, util::json::Dict{{"statistics", totalScore},
                                      {"component-statistics", jsonScores}});
      for (auto res : resultGraphs) {
        gout.printLatLng(*res, &out);
      }
      out.flush();
    } else {
      util::geo::output::GeoJsonOutput out(std::cout);

      for (auto res : resultGraphs) {
        gout.printLatLng(*res, &out);
      }
      out.flush();
    }
  }

  return 0;
}
