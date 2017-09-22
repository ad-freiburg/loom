// Copyright 2017, University of Freiburg,
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
  config.add_options()(
      "grid-cell-width",
      opts::value<double>(&(cfg->gridCellWidth))->default_value(50.0),
      "grid cell width")
    (
      "print-mode",
      opts::value<std::string>(&(cfg->printMode))->default_value("transitgraph"),
      "print mode: transitgraph, gridgraph")
        (
      "vert-pen",
      opts::value<double>(&(cfg->vertPen))->default_value(3),
      "penalty for vertical edges")(
      "hori-pen",
      opts::value<double>(&(cfg->horiPen))->default_value(3),
      "penalty for horicontal edges")(
      "diag-pen",
      opts::value<double>(&(cfg->diagPen))->default_value(5),
      "penalty for diagonal edges");

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
