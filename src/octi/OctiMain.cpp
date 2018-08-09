// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include "octi/Octilinearizer.h"
#include "octi/config/ConfigReader.h"
#include "octi/combgraph/CombGraph.h"
#include "octi/transitgraph/TransitGraph.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/Geo.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

using octi::Octilinearizer;

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

  T_START(octilinearize);
  Octilinearizer oct;
  TransitGraph res = oct.draw(&tg, &gg, cfg.pens);
  std::cerr << " octilinearized input in " << T_STOP(octilinearize) << "ms"
            << std::endl;

  if (cfg.printMode == "gridgraph") {
    out.print(*gg, std::cout);
  } else {
    out.print(res, std::cout);
  }

  return (0);
}
