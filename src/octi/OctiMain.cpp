// Copyright 2017
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "octi/config/ConfigReader.h"
#include "octi/graph/Graph.h"
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/graph/Dijkstra.h"

using std::string;
using namespace octi;

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
   *   gridgraph::GridGraph g(util::geo::Box(util::geo::Point(0, 0),
   * util::geo::Point(100, 100)), 10);
   *   g.balance();
   *
   *   util::graph::Dijkstra dijkstra;
   *
   *   std::list<util::graph::Edge<gridgraph::NodePL, gridgraph::EdgePL>*> res;
   *
   *   dijkstra.shortestPath(g.getNode(0, 0), g.getNode(5, 8), &res);
   *
   *   for (auto e : res) {
   *     e->pl().addRoute("A");
   *   }
   *
   *   res.clear();
   *   g.balance();
   *
   *   dijkstra.shortestPath(g.getNode(0, 8), g.getNode(9, 0), &res);
   *
   *   for (auto e : res) {
   *     e->pl().addRoute("B");
   *   }
   *
   *   res.clear();
   *   g.balance();
   *
   *   out.print(g);
   */

  graph::Graph tg(&(std::cin));

  std::map<util::graph::Node<octi::graph::NodePL, octi::graph::EdgePL>,
           util::graph::Node<octi::gridgraph::NodePL, octi::gridgraph::EdgePL>>
      m;

  std::set<util::graph::Node<octi::gridgraph::NodePL, octi::gridgraph::EdgePL>> used;

  gridgraph::GridGraph g(tg.getBBox(), 500);

  for (auto n : *tg.getNodes()) {
    g.get    
  }

  g.balance();

  out.print(g);

  return (0);
}
