// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <chrono>
#include "octi/Octilinearizer.h"
#include "octi/config/ConfigReader.h"
#include "octi/graph/CombEdgePL.h"
#include "octi/graph/CombNodePL.h"
#include "octi/graph/Graph.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

using octi::gridgraph::GridGraph;
using util::graph::Node;
using octi::graph::NodePL;
using octi::graph::EdgePL;
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

  /*
   *   {
   *     gridgraph::GridGraph g(
   *         util::geo::Box(util::geo::Point(0, 0), util::geo::Point(100, 100)),
   * 10);
   *
   *     util::graph::Dijkstra dijkstra;
   *     std::list<GridEdge*> res;
   *
   *     dijkstra.shortestPath(g.getNode(4, 2), g.getNode(3, 8), &res);
   *
   *     for (auto f : res) {
   *       f->pl().addResidentEdge(0);
   *     }
   *
   *     for (auto f : res) {
   *       g.balanceEdge(f->getFrom()->pl().getParent(),
   *                     f->getTo()->pl().getParent());
   *     }
   *
   *     res.clear();
   *     dijkstra.shortestPath(g.getNode(8, 2), g.getNode(3, 8), &res);
   *
   *     for (auto f : res) {
   *       f->pl().addResidentEdge(0);
   *     }
   *
   *     for (auto f : res) {
   *       g.balanceEdge(f->getFrom()->pl().getParent(),
   *                     f->getTo()->pl().getParent());
   *     }
   *
   *     out.print(g);
   *     exit(0);
   *   }
   */
  std::chrono::steady_clock::time_point begin;
  std::chrono::steady_clock::time_point end;

  std::cerr << "Reading graph file... ";
  begin = std::chrono::steady_clock::now();
  TransitGraph tg;
  GridGraph* gg;
  if (cfg.fromDot) tg.readFromDot(&(std::cin));
  else tg.readFromJson(&(std::cin));
  end = std::chrono::steady_clock::now();
  std::cerr << " done (" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms)"<< std::endl;


  std::cerr << "Planarize graph... ";
  begin = std::chrono::steady_clock::now();
  tg.topologizeIsects();
  end = std::chrono::steady_clock::now();
  std::cerr << " done (" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms)"<< std::endl;

  auto totalBegin = std::chrono::steady_clock::now();
  Octilinearizer oct;
  TransitGraph res = oct.draw(&tg, &gg, cfg.pens);
  auto totalend = std::chrono::steady_clock::now();
  std::cerr << " octilinearized input in " << std::chrono::duration_cast<std::chrono::milliseconds>(totalend - totalBegin).count() << "ms"<< std::endl;

  if (cfg.printMode == "gridgraph") {
    out.print(*gg, std::cout);
  } else {
    out.print(res, std::cout);
  }

  return (0);
}
