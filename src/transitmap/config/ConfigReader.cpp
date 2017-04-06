// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <boost/program_options.hpp>
#include <iostream>
#include <float.h>
#include <string>
#include <exception>
#include "pbutil/log/Log.h"
#include "./ConfigReader.h"

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
     opts::bool_switch(&(cfg->renderStations))
      ->default_value(true),
      "render station geometries")
    ("render-station-names",
     opts::bool_switch(&(cfg->renderStationNames))
      ->default_value(false),
      "render station name")
    ("render-node-fronts",
     opts::bool_switch(&(cfg->renderNodeFronts))
      ->default_value(false),
      "render node fronts as red lines")
    ("no-optim,N",
      opts::bool_switch(&(cfg->noOptim))
      ->default_value(false),
      "disable line-ordering optimization")
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
    ("input-smoothing",
      opts::value<double>(&(cfg->inputSmoothing))
      ->default_value(3),
      "level of input-data smoothing")
    ("bezier-prec",
      opts::value<double>(&(cfg->innerGeometryPrecision))
      ->default_value(3),
      "rendering precision of inner node geometries")
    ("in-station-cross-penalty-factor",
      opts::value<double>(&(cfg->inStationCrossPenalty))
      ->default_value(3),
      "penalty factor during optimization for crossings that occur in stations")
    ("glpk-time-limit",
      opts::value<int>(&(cfg->glpkTimeLimit))
      ->default_value(60000 + 60000),
      "GLPK: overall time limit for search, in ms")
    ("glpk-proximity-search-time-limit",
      opts::value<int>(&(cfg->glpkPSTimeLimit))
      ->default_value(60000),
      "GLPK: time limit for proximit search heuristig")
    ("glpk-use-proximity-search",
      opts::bool_switch(&(cfg->useGlpkProximSearch))
      ->default_value(true),
      "GLPK: use proximity search heuristic")
    ("glpk-use-feasibility-pump",
      opts::bool_switch(&(cfg->useGlpkProximSearch))
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

