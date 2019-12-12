// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include "octi/Octilinearizer.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/config/ConfigReader.h"
#include "octi/gridgraph/GridGraph.h"
#include "shared/linegraph/LineGraph.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

using octi::Octilinearizer;
using util::geo::dist;

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
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  util::geo::output::GeoGraphJsonOutput out;

  std::cerr << "Reading graph file... ";
  T_START(read);
  LineGraph tg;
  GridGraph* gg;
  if (cfg.fromDot)
    tg.readFromDot(&(std::cin));
  else
    tg.readFromJson(&(std::cin));
  std::cerr << " done (" << T_STOP(read) << "ms)" << std::endl;

  std::cerr << "Planarize graph... ";
  T_START(planarize);
  tg.topologizeIsects();
  std::cerr << " done (" << T_STOP(planarize) << "ms)" << std::endl;

  double avgDist = avgStatDist(tg);
  std::cerr << "Average adj. node distance is " << avgDist << std::endl;

  Octilinearizer oct;
  LineGraph res;

  double gridSize;

  if (util::trim(cfg.gridSize).back() == '%') {
    double perc = atof(cfg.gridSize.c_str()) / 100;
    gridSize = avgDist * perc;
    std::cerr << "Grid size " << gridSize << " (" << perc * 100 << "%)\n";
  } else {
    gridSize = atof(cfg.gridSize.c_str());
    std::cerr << "Grid size " << gridSize << "\n";
  }

  if (cfg.optMode == "ilp") {
    T_START(octi);
    double sc =
        oct.drawILP(&tg, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                    cfg.deg2Heur, cfg.maxGrDist, cfg.ilpNoSolve, cfg.ilpPath);
    std::cerr << "Octilinearized using ILP in " << T_STOP(octi) << " ms, score "
              << sc << std::endl;
  } else if ((cfg.optMode == "heur")) {
    T_START(octi);
    double sc = oct.draw(&tg, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                         cfg.deg2Heur, cfg.maxGrDist, cfg.restrLocSearch,
                         cfg.enfGeoPen);
    std::cerr << "Octilinearized using heur approach in " << T_STOP(octi)
              << " ms, score " << sc << std::endl;
  }

  if (cfg.printMode == "gridgraph") {
    out.print(*gg, std::cout);
  } else {
    out.print(res, std::cout);
  }

  return (0);
}
