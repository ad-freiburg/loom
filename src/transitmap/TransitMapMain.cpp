// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "shared/rendergraph/Penalties.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/ConfigReader.cpp"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/GraphBuilder.h"
#include "transitmap/output/SvgRenderer.h"
#include "util/log/Log.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  transitmapper::config::Config cfg;

  transitmapper::config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  LOGTO(DEBUG, std::cerr) << "Reading graph...";
  shared::rendergraph::RenderGraph g(cfg.lineWidth, cfg.lineSpacing);
  transitmapper::graph::GraphBuilder b(&cfg);

  if (cfg.fromDot) {
    g.readFromDot(&std::cin, cfg.inputSmoothing);
  } else {
    g.readFromJson(&std::cin, cfg.inputSmoothing);
  }

  g.smooth();

  b.writeNodeFronts(&g);

  b.expandOverlappinFronts(&g);

  // find expanded node fronts that form a node and replace them with a
  // single node
  g.createMetaNodes();

  if (cfg.renderMethod == "svg") {
    LOGTO(DEBUG, std::cerr) << "Outputting to SVG ...";
    transitmapper::output::SvgRenderer svgOut(&std::cout, &cfg);
    svgOut.print(g);
  } else {
    LOG(ERROR) << "Unknown render method " << cfg.renderMethod;
    exit(1);
  }

  return (0);
}
