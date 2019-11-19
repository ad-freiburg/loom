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
#include "shared/transitgraph/TransitGraph.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

using octi::Octilinearizer;
using util::geo::dist;

// _____________________________________________________________________________
double avgStatDist(const TransitGraph& g) {
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
  TransitGraph tg;
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
  TransitGraph res;

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
    double sc = oct.drawILP(&tg, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                            cfg.deg2Heur);
    std::cerr << "Octilinearized using ILP in " << T_STOP(octi) << " ms, score "
              << sc << std::endl;
  } else if ((cfg.optMode == "heur")) {
    T_START(octi);
    double sc = oct.draw(&tg, &res, &gg, cfg.pens, gridSize, cfg.borderRad,
                         cfg.deg2Heur);
    std::cerr << "Octilinearized using heur approach in " << T_STOP(octi)
              << " ms, score " << sc << std::endl;
  } else if ((cfg.optMode == "eval")) {
    TransitGraph resIlp, resHeur;
    T_START(octi_ilp);
    double ilpSc = oct.drawILP(&tg, &resIlp, &gg, cfg.pens, gridSize,
                               cfg.borderRad, cfg.deg2Heur);
    auto ilpT = T_STOP(octi_ilp);
    T_START(octi_heur);
    double heurSc = oct.draw(&tg, &resHeur, &gg, cfg.pens, gridSize,
                             cfg.borderRad, cfg.deg2Heur);
    auto heurT = T_STOP(octi_heur);

    std::cerr << "\nOctilinearized using eval approach: " << std::endl;
    std::cerr << "  ILP  target value: " << ilpSc << std::endl;
    std::cerr << "  HEUR target value: " << heurSc << " (+"
              << std::setprecision(2) << std::fixed
              << (((heurSc - ilpSc) / ilpSc) * 100) << "%)" << std::endl;
    std::cerr << "  ILP  solve time: " << std::setprecision(2) << std::fixed
              << ilpT << " ms" << std::endl;
    std::cerr << "  HEUR solve time: " << std::setprecision(2) << std::fixed
              << heurT << " ms" << std::endl;

    std::ofstream of;
    of.open(cfg.evalPath + "/res_ilp" + cfg.evalSuff + ".json");
    out.print(resIlp, of);
    of.flush();
    of.close();
    of.open(cfg.evalPath + "/res_heur" + cfg.evalSuff + ".json");
    out.print(resHeur, of);
    of.flush();
    of.close();
  }

  if ((cfg.optMode != "eval")) {
    if (cfg.printMode == "gridgraph") {
      out.print(*gg, std::cout);
    } else {
      out.print(res, std::cout);
    }
  }

  return (0);
}
