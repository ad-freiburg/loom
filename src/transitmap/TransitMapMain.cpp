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
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace transitmapper;
using namespace util::geo;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  config::Config cfg;

  config::ConfigReader cr;
  cr.read(&cfg, argc, argv);

  LOG(DEBUG) << "Reading graph...";
  shared::rendergraph::RenderGraph g(cfg.lineWidth, cfg.lineSpacing);
  transitmapper::graph::GraphBuilder b(&cfg);

  g.readFromJson(&std::cin, cfg.inputSmoothing);

  LOG(DEBUG) << "Creating node fronts...";
  b.writeMainDirs(&g);

  if (cfg.expandFronts) b.expandOverlappinFronts(&g);

  g.createMetaNodes();

  if (cfg.renderMethod == "svg") {
    std::string path = cfg.outputPath;
    LOG(DEBUG) << "Outputting to SVG " << path << " ...";
    std::ofstream o;
    o.open(path);
    output::SvgRenderer svgOut(&o, &cfg);
    svgOut.print(g);
  }

  if (cfg.renderMethod == "svg_sep") {
    {
      std::string path = cfg.outputPath + "/edges.svg";
      cfg.renderEdges = 1;
      cfg.renderStations = 0;
      cfg.renderLabels = 0;
      cfg.renderNodeConnections = 0;
      LOG(DEBUG) << "Outputting edge SVG to " << path << " ...";
      std::ofstream o;
      o.open(path);
      output::SvgRenderer svgOut(&o, &cfg);
      svgOut.print(g);
    }

    {
      std::string path = cfg.outputPath + "/nodes.svg";
      cfg.renderEdges = 0;
      cfg.renderStations = 0;
      cfg.renderLabels = 0;
      cfg.renderNodeConnections = 1;
      LOG(DEBUG) << "Outputting node SVG to " << path << " ...";
      std::ofstream o;
      o.open(path);
      output::SvgRenderer svgOut(&o, &cfg);
      svgOut.print(g);
    }

    {
      std::string path = cfg.outputPath + "/stations.svg";
      cfg.renderEdges = 0;
      cfg.renderStations = 1;
      cfg.renderLabels = 0;
      cfg.renderNodeConnections = 0;
      LOG(DEBUG) << "Outputting station SVG to " << path << " ...";
      std::ofstream o;
      o.open(path);
      output::SvgRenderer svgOut(&o, &cfg);
      svgOut.print(g);
    }

    {
      std::string path = cfg.outputPath + "/labels.svg";
      cfg.renderEdges = 0;
      cfg.renderStations = 0;
      cfg.renderLabels = 1;
      cfg.renderNodeConnections = 0;
      LOG(DEBUG) << "Outputting label SVG to " << path << " ...";
      std::ofstream o;
      o.open(path);
      output::SvgRenderer svgOut(&o, &cfg);
      svgOut.print(g);
    }
  }

  return (0);
}
