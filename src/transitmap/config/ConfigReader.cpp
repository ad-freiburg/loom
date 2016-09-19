// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <boost/program_options.hpp>
#include <iostream>
#include <float.h>
#include <string>
#include <exception>
#include "./../../log/Log.h"
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
    ("verbose,v", "verbosity level");

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
  ;

  opts::options_description cmdline_options;
  cmdline_options.add(config).add(generic);

  opts::options_description config_file_options;
  config_file_options.add(config);

  opts::options_description visibleDesc("Allowed options");
  visibleDesc.add(generic).add(config);
  opts::variables_map vm;

  try {
    opts::store(opts::command_line_parser(argc, argv).options(cmdline_options)
              .run(), vm);
    opts::notify(vm);
  } catch (exception e) {
    LOG(ERROR) << e.what();
    std::cout << visibleDesc << "\n";
    exit(1);
  }

  if (vm.count("help")) {
    std::cout << "transitmapper\n"
      << VERSION_STR << "\n\n"
      << visibleDesc << "\n";
    exit(0);
  }

  if (vm.count("version")) {
    std::cout << "\n" << VERSION_STR << "\n";
    exit(0);
  }
}

