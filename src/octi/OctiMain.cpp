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
#include "octi/gridgraph/GridGraph.h"
#include "util/geo/Geo.h"
#include "util/graph/Dijkstra.h"
#include "util/geo/output/GeoJsonOutput.h"

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

  gridgraph::GridGraph g(util::geo::Box(util::geo::Point(0, 0), util::geo::Point(100, 100)), 10);
  util::geo::output::GeoJsonOutput out;

  util::graph::Dijkstra dijkstra;

  std::list<util::graph::Edge<gridgraph::NodePL, gridgraph::EdgePL>*> res;

  dijkstra.shortestPath(g.getNode(0, 0), g.getNode(5, 8), &res);

  for (auto e : res) {
    e->pl().addRoute("A");
  }

  res.clear();

  dijkstra.shortestPath(g.getNode(0, 8), g.getNode(9, 0), &res);

  for (auto e : res) {
    e->pl().addRoute("B");
  }


  out.print(g);

  return (0);
}
