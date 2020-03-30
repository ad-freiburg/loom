// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <boost/program_options.hpp>
#include <exception>
#include <iostream>
#include <string>
#include "octi/config/ConfigReader.h"
#include "util/log/Log.h"

using octi::config::ConfigReader;
namespace opts = boost::program_options;

using std::string;
using std::exception;
using std::vector;

// _____________________________________________________________________________
ConfigReader::ConfigReader() {}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) const {
  string VERSION_STR = " - unversioned - ";

  opts::options_description generic("General");
  generic.add_options()("version", "output version")(
      "help,?", "show this message")("verbose,v", "verbosity level");

  opts::options_description config("Output");
  config.add_options()
    (
      "optim-mode,o",
      opts::value<std::string>(&(cfg->optMode))->default_value("heur"),
      "optimization mode, either 'heur' (fast) or 'ilp' (very slow)")
    (
      "ilp-no-solve",
      opts::bool_switch(&(cfg->ilpNoSolve))->default_value(false),
      "if set, the ILP is not solved, only written to file")
    ("ilp-time-limit",
      opts::value<int>(&(cfg->ilpTimeLimit))
      ->default_value(60),
      "ILP solver time limit (in seconds), negative value means infinity")
    ("ilp-solver",
      opts::value<std::string>(&(cfg->ilpSolver))
      ->default_value("gurobi"),
      "The preferred solver library to use, will fall back if library is not available.")
    (
      "ilp-out",
      opts::value<std::string>(&(cfg->ilpPath))->default_value("prob.mps"),
      "path to output the ILP to, a first feasible solution will be written "
      "to <basename>.mst.")
    (
      "obstacles",
      opts::value<std::string>(&(cfg->obstaclePath))->default_value(""),
      "GeoJSON file containing obstacle polygons")
    (
      "from-dot,D",
      opts::bool_switch(&(cfg->fromDot))->default_value(false),
      "input is in dot format")
    (
      "deg2-heur",
      opts::bool_switch(&(cfg->deg2Heur))->default_value(false),
      "contract degree 2 nodes and re-insert them equidistantly")
    (
      "enf-geo-pen",
      opts::value<double>(&(cfg->enfGeoPen))->default_value(0),
      "enforce lines to roughly follow original geographical course")
    (
      "max-grid-dist",
      opts::value<double>(&(cfg->maxGrDist))->default_value(3),
      "max grid distance radius for station candidates")
    (
      "restr-loc-search",
      opts::bool_switch(&(cfg->restrLocSearch))->default_value(false),
      "restrict local search to max grid distance")
    (
      "density-pen",
      opts::value<double>(&(cfg->pens.densityPen))->default_value(0),
      "penalty factor for re-inserted contracted stations that are too near, a reasonable value is e.g. 5. Only works with optim mode 'heur'!")
    (
      "grid-size,g",
      opts::value<std::string>(&(cfg->gridSize))->default_value("100%"),
      "grid cell length, either as exact value (like '500') or as percentage of average adj. station distance (like '75%')")
    (
      "border-rad",
      opts::value<double>(&(cfg->borderRad))->default_value(45),
      "border rad")
    (
      "print-mode",
      opts::value<std::string>(&(cfg->printMode))->default_value("transitgraph"),
      "print mode: transitgraph, gridgraph")
    (
      "vert-pen",
      opts::value<double>(&(cfg->pens.verticalPen))->default_value(0),
      "penalty for vertical edges")
    (
      "hori-pen",
      opts::value<double>(&(cfg->pens.horizontalPen))->default_value(0),
      "penalty for horicontal edges")
    (
      "diag-pen",
      opts::value<double>(&(cfg->pens.diagonalPen))->default_value(.5),
      "penalty for diagonal edges")
    (
      "pen-180",
      opts::value<double>(&(cfg->pens.p_0))->default_value(0),
      "penalty for 180 deg traversal")
    (
      "pen-135",
      opts::value<double>(&(cfg->pens.p_135))->default_value(1),
      "penalty for 135 deg traversal")
    (
      "pen-90",
      opts::value<double>(&(cfg->pens.p_90))->default_value(1.5),
      "penalty for 90 deg traversal")
    (
      "pen-45",
      opts::value<double>(&(cfg->pens.p_45))->default_value(2),
      "penalty for 45 deg traversal")
    ;

  opts::options_description cmdlineOptions;
  cmdlineOptions.add(config).add(generic);

  opts::options_description visibleDesc("Allowed options");
  visibleDesc.add(generic).add(config);
  opts::variables_map vm;

  try {
    opts::store(
        opts::command_line_parser(argc, argv).options(cmdlineOptions).run(),
        vm);
    opts::notify(vm);
  } catch (exception e) {
    LOG(ERROR) << e.what() << std::endl;
    std::cout << visibleDesc << "\n";
    exit(1);
  }

  if (vm.count("help")) {
    std::cout << argv[0] << " [options] <input-feed>\n"
              << VERSION_STR << "\n\n"
              << visibleDesc << "\n";
    exit(0);
  }

  if (vm.count("version")) {
    std::cout << "\n" << VERSION_STR << "\n";
    exit(0);
  }
}
