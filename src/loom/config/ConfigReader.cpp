// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <boost/program_options.hpp>
#include <iostream>
#include <float.h>
#include <string>
#include <exception>
#include "util/log/Log.h"
#include "loom/config/ConfigReader.h"
#include "loom/config/LoomConfig.h"

using loom::config::ConfigReader;
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
      "output path");

  opts::options_description config("Output");
  config.add_options()
    ("name",
      opts::value<std::string>(&(cfg->name))
      ->default_value("shinygraph"),
      "name of transit graph")
    ("untangle",
     opts::value<bool>(&(cfg->untangleGraph))
      ->default_value(true),
      "untangle line graph")
    ("output-stats",
      opts::value<bool>(&(cfg->outputStats))
      ->default_value(false),
      "print some more stats")
    ("prune",
      opts::value<bool>(&(cfg->pruneGraph))
      ->default_value(true),
      "prune and cut the input graph prior to optimization")
    ("optim-method",
      opts::value<std::string>(&(cfg->optimMethod))
      ->default_value("comb"),
      "optimization method to use, possible values: ilp_impr, ilp, hillc, comb")
    ("same-seg-cross-penalty-factor",
      opts::value<double>(&(cfg->crossPenMultiSameSeg))
      ->default_value(4),
      "penalty factor during optimization for crossings that occur between two lines that travel 2 same segments")
    ("separation-penalty-factor",
      opts::value<double>(&(cfg->separationPenWeight))
      ->default_value(3),
      "penalty factor during optimization for separations")
    ("in-station-crossing-penalty-factor-same-seg",
      opts::value<double>(&(cfg->stationCrossWeightSameSeg))
      ->default_value(12),
      "penalty factor during optimization for crossings in station with degree > 2 that occur between two lines that travel 2 same segments, only applies if degree-based penalty is enabled")
   ("in-station-crossing-penalty-factor-diff-seg",
      opts::value<double>(&(cfg->stationCrossWeightDiffSeg))
      ->default_value(3),
      "penalty factor during optimization for crossings in station with degree > 2 that occur between two lines that travel 2 same segments, only applies if degree-based penalty is enabled")
    ("in-station-separation-penalty-factor",
      opts::value<double>(&(cfg->stationSeparationWeight))
      ->default_value(9),
      "penalty factor during optimization for separations in station with degree > 2, only applies if degree-based penalty is enabled")
    ("diff-seg-cross-penalty-factor",
      opts::value<double>(&(cfg->crossPenMultiDiffSeg))
      ->default_value(1),
      "penalty factor during optimization for crossings that occur between two lines that travel through 3 segments")
    ("ilp-time-limit",
      opts::value<int>(&(cfg->ilpTimeLimit))
      ->default_value(-1),
      "ILP solver time limit (in seconds), negative value means infinity")
    ("ilp-solver",
      opts::value<std::string>(&(cfg->ilpSolver))
      ->default_value("gurobi"),
      "The preferred solver library to use, will fall back if library is not available.")
    ("mps-output-path",
      opts::value<std::string>(&(cfg->MPSOutputPath))
      ->default_value(""),
      "output path for ILP, printed in MPS format. If empty, no output will be generated.")
    ("optim-runs",
      opts::value<size_t>(&(cfg->optimRuns))
      ->default_value(1),
      "number of optimization runs (use for evaluation, avg stats will be printed if > 1)")
    ("dbg-output-path",
      opts::value<std::string>(&(cfg->dbgPath))
      ->default_value("."),
      "folder the debug information will be written to")
    ("output-optgraph",
      opts::value<bool>(&(cfg->outOptGraph))
      ->default_value(false),
      "write optimization graph to debug path (--dbg-output-path)")
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

