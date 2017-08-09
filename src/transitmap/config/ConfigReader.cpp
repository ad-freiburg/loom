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
    ("render-stations",
     opts::value<bool>(&(cfg->renderStations))
      ->default_value(true),
      "render station geometries")
    ("dont-expand-stations",
     opts::value<bool>(&(cfg->dontExpandStations))
      ->default_value(true),
      "dont expand node fronts for station rendering")
    ("render-station-names",
     opts::value<bool>(&(cfg->renderStationNames))
      ->default_value(false),
      "render station name")
    ("render-dir-markers",
     opts::value<bool>(&(cfg->renderDirMarkers))
      ->default_value(false),
      "render direction markers (experimental!)")
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
    ("create-core-optim-graph",
     opts::value<bool>(&(cfg->createCoreOptimGraph))
      ->default_value(true),
      "simplify optimization graph prior to optimization")
    ("use-cbc",
     opts::value<bool>(&(cfg->useCbc))
      ->default_value(true),
      "use the COIN-CBC solver for optimization instead of GLPK")
    ("render-node-circles",
     opts::value<bool>(&(cfg->renderNodeCircles))
      ->default_value(false),
      "mark node areas with background grey")
    ("render-node-polygons",
     opts::value<bool>(&(cfg->renderNodePolygons))
      ->default_value(false),
      "render node polygons")
    ("no-optim,N",
      opts::value<bool>(&(cfg->noOptim))
      ->default_value(false),
      "disable line-ordering optimization")
    ("splitting-optim",
      opts::value<bool>(&(cfg->splittingOpt))
      ->default_value(false),
      "enable splitting optimization")
    ("output-stats",
      opts::value<bool>(&(cfg->outputStats))
      ->default_value(false),
      "print some more graph stats to stdout")
    ("render-stats",
      opts::value<bool>(&(cfg->renderStats))
      ->default_value(false),
      "render stats to output")
    ("collapse-line-partners",
      opts::value<bool>(&(cfg->collapseLinePartners))
      ->default_value(true),
      "collapse line partners")
    ("simple-station-render-heur",
      opts::value<bool>(&(cfg->simpleRenderForTwoEdgeNodes))
      ->default_value(true),
      "if a station only has 2 incident edges, dont calculate hull")
    ("optim-method",
      opts::value<std::string>(&(cfg->optimMethod))
      ->default_value("ilp_impr"),
      "optimization method to use, possible values: ilp_impr, ilp, hillc")
    ("projection",
      opts::value<std::string>(&(cfg->projectionString))
      ->default_value("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"),
      "map projection as proj4 string")
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
    ("same-seg-cross-penalty-factor",
      opts::value<double>(&(cfg->crossPenMultiSameSeg))
      ->default_value(4),
      "penalty factor during optimization for crossings that occur between two lines that travel 2 same segments")
    ("splitting-penalty-factor",
      opts::value<double>(&(cfg->splitPenWeight))
      ->default_value(3),
      "penalty factor during optimization for splittings")
    ("in-station-crossing-penalty-factor",
      opts::value<double>(&(cfg->stationCrossWeight))
      ->default_value(3),
      "penalty factor during optimization for splittings")
    ("diff-seg-cross-penalty-factor",
      opts::value<double>(&(cfg->crossPenMultiDiffSeg))
      ->default_value(1),
      "penalty factor during optimization for crossings that occur between two lines that travel through 3 segments")
    ("glpk-time-limit",
      opts::value<int>(&(cfg->glpkTimeLimit))
      ->default_value(60000 + 60000),
      "GLPK: overall time limit for search, in ms")
    ("glpk-proximity-search-time-limit",
      opts::value<int>(&(cfg->glpkPSTimeLimit))
      ->default_value(60000),
      "GLPK: time limit for proximit search heuristig")
    ("glpk-use-proximity-search",
      opts::value<bool>(&(cfg->useGlpkProximSearch))
      ->default_value(true),
      "GLPK: use proximity search heuristic")
    ("glpk-use-feasibility-pump",
      opts::value<bool>(&(cfg->useGlpkProximSearch))
      ->default_value(true),
      "GLPK: use feasibility pump heuristic")
    ("glpk-mps-output-path",
      opts::value<std::string>(&(cfg->glpkMPSOutputPath))
      ->default_value(""),
      "output path for ILP, printed in MPS format. If empty, no output will be generated.")
    ("glpk-h-output-path",
      opts::value<std::string>(&(cfg->glpkHOutputPath))
      ->default_value(""),
      "output path for ILP, printed in \"human readable\" format. If empty, no output will be generated.")
    ("glpk-solution-output-path",
      opts::value<std::string>(&(cfg->glpkSolutionOutputPath))
      ->default_value(""),
      "output path for ILP solution. If empty, no output will be generated.")
    ("coin-cbc-binary",
      opts::value<std::string>(&(cfg->coinCbcBinary))
     ->default_value("/home/patrick/repos/Cbc-2.9/build/bin/cbc"),
      "coin cbc binary path")
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

