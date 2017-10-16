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
#include "util/geo/output/GeoJsonOutput.h"
#include "util/log/Log.h"

using std::string;
using namespace octi;

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

  util::geo::output::GeoJsonOutput out;

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
  auto totalBegin = std::chrono::steady_clock::now();

  std::cerr << "Reading graph file... ";
  begin = std::chrono::steady_clock::now();
  TransitGraph tg(&(std::cin));
  end = std::chrono::steady_clock::now();
  std::cerr << " done (" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms)"<< std::endl;

  std::cerr << "Planarize graph... ";
  begin = std::chrono::steady_clock::now();
  tg.topologizeIsects();
  end = std::chrono::steady_clock::now();
  std::cerr << " done (" << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count() << "ms)"<< std::endl;

  Octilinearizer oct;
  TransitGraph res = oct.draw(&tg, cfg.pens);

  if (cfg.printMode == "gridgraph") {
    /*
     * auto nodes = *g.getNodes();
     * for (auto n : nodes) {
     *   if (n->getAdjList().size() == 16) {
     *     g.deleteNode(n);
     *   } else {
     *     for (auto e : n->getAdjListOut()) {
     *       if (n->pl().getParent() == e->getOtherNode(n)->pl().getParent()) {
     *         g.deleteEdge(n, e->getOtherNode(n));
     *       }
     *     }
     *   }
     * }
     */
    //out.print(g);
    //exit(0);
  } else {
    out.print(res);
  }

  auto totalend = std::chrono::steady_clock::now();
  std::cerr << " octilinearized input in " << std::chrono::duration_cast<std::chrono::milliseconds>(totalend - totalBegin).count() << "ms"<< std::endl;

  return (0);
}
