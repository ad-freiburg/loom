// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <boost/program_options.hpp>
#include <iostream>
#include <float.h>
#include <string>
#include <exception>
#include "util/log/Log.h"
#include "transitmap/config/ConfigReader.h"

using transitmapper::config::ConfigReader;
namespace opts = boost::program_options;
using std::string;
using std::exception;
using std::vector;

// _____________________________________________________________________________
ConfigReader::ConfigReader() {
}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) const {
  string VERSION_STR = " - unversioned - ";

  opts::options_description generic("General");
  generic.add_options()
    ("version", "output version")
    ("help,?", "show this message")
    ("verbose,v", "verbosity level")
    ("output,o",
      opts::value<std::string>(&(cfg->outputPath))
      ->default_value("./transitmap.svg"),
      "output path")
    ("render-engine",
      opts::value<std::string>(&(cfg->renderMethod))
      ->default_value("svg"),
      "output method, atm only \"svg\" is supported");

  opts::options_description config("Output");
  config.add_options()
    ("name",
      opts::value<std::string>(&(cfg->name))
      ->default_value("shinygraph"),
      "name of transit graph")
    ("line-width",
      opts::value<double>(&(cfg->lineWidth))
      ->default_value(20),
      "width of a single rendered line")
    ("line-spacing",
      opts::value<double>(&(cfg->lineSpacing))
      ->default_value(10),
      "spacing between two rendered lines")
    ("outline-width",
      opts::value<double>(&(cfg->outlineWidth))
      ->default_value(2),
      "default width of line outline")
    ("line-label-textsize",
      opts::value<double>(&(cfg->lineLabelSize))
      ->default_value(40),
      "text size for line labels")
    ("station-label-textsize",
      opts::value<double>(&(cfg->stationLabelSize))
      ->default_value(60),
      "text size for station labels")
    ("render-stations",
     opts::value<bool>(&(cfg->renderStations))
      ->default_value(true),
      "render station geometries")
    ("render-labels",
     opts::value<bool>(&(cfg->renderLabels))
      ->default_value(false),
      "render labels")
    ("tight-stations",
     opts::value<bool>(&(cfg->tightStations))
      ->default_value(false),
      "dont expand station nodes at all")
    ("render-dir-markers",
     opts::value<bool>(&(cfg->renderDirMarkers))
      ->default_value(false),
      "render direction markers")
    ("render-node-connections",
     opts::value<bool>(&(cfg->renderNodeConnections))
      ->default_value(true),
      "render inner node connections")
    ("render-edges",
     opts::value<bool>(&(cfg->renderEdges))
      ->default_value(true),
      "render edges")
    ("expand-fronts",
     opts::value<bool>(&(cfg->expandFronts))
      ->default_value(true),
      "expand fronts to make room for inner geometries")
    ("render-node-fronts",
     opts::value<bool>(&(cfg->renderNodeFronts))
      ->default_value(false),
      "render node fronts as red lines")
    ("render-node-circles",
     opts::value<bool>(&(cfg->renderNodeCircles))
      ->default_value(false),
      "mark node areas with background grey")
    ("render-node-polygons",
     opts::value<bool>(&(cfg->renderNodePolygons))
      ->default_value(false),
      "render node polygons")
    ("render-stats",
      opts::value<bool>(&(cfg->renderStats))
      ->default_value(false),
      "render stats to output")
    ("simple-station-render-heur",
      opts::value<bool>(&(cfg->simpleRenderForTwoEdgeNodes))
      ->default_value(true),
      "if a station only has 2 incident edges, dont calculate hull")
    ("resolution",
      opts::value<double>(&(cfg->outputResolution))
      ->default_value(0.1),
      "output resolution (output pixel size / projection pixel size)")
    ("padding",
      opts::value<double>(&(cfg->outputPadding))
      ->default_value(-1),
      "padding of output. If < 0, automatic.")
    ("input-smoothing",
      opts::value<double>(&(cfg->inputSmoothing))
      ->default_value(3),
      "level of input-data smoothing")
    ("bezier-prec",
      opts::value<double>(&(cfg->innerGeometryPrecision))
      ->default_value(3),
      "rendering precision of inner node geometries")
    ("world-file-path",
      opts::value<std::string>(&(cfg->worldFilePath))
      ->default_value(""),
      "if set, output world file of rendered map to this location")
  ;

  opts::options_description cmdlineOptions;
  cmdlineOptions.add(config).add(generic);

  opts::options_description visibleDesc("Allowed options");
  visibleDesc.add(generic).add(config);
  opts::variables_map vm;

  try {
    opts::store(opts::command_line_parser(argc, argv).options(cmdlineOptions)
        .run(), vm);
    opts::notify(vm);

    if (cfg->outputPadding < 0) {
      cfg->outputPadding = (cfg->lineWidth + cfg->lineSpacing);
    }
  } catch (exception e) {
    LOG(ERROR) << e.what() << std::endl;
    std::cout << visibleDesc << "\n";
    exit(1);
  }

  if (vm.count("help")) {
    std::cout << argv[0] << " [options]\n"
      << VERSION_STR << "\n\n"
      << visibleDesc << "\n";
    exit(0);
  }

  if (vm.count("version")) {
    std::cout << "\n" << VERSION_STR << "\n";
    exit(0);
  }
}

