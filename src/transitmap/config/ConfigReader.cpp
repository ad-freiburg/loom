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
    ("ordering-optim,O",
      opts::value<size_t>(&(cfg->optimIterations))
      ->default_value(10),
      "number of edge order optimization iterations")
    ("projection",
      opts::value<std::string>(&(cfg->projectionString))
      ->default_value("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"),
      "map projection as proj4 string")
    ("resolution",
      opts::value<double>(&(cfg->outputResolution))
      ->default_value(0.1),
      "output resolution (output pixel size / projection pixel size)")
  ;

  opts::options_description positional("Positional arguments");
  positional.add_options()
    ("input-feed",
      opts::value<std::string>(&(cfg->inputFeedPath)),
      "path to an (unzipped) GTFS feed");

  opts::positional_options_description positionalOptions;
  positionalOptions.add("input-feed", 1);

  opts::options_description cmdlineOptions;
  cmdlineOptions.add(config).add(generic).add(positional);

  opts::options_description visibleDesc("Allowed options");
  visibleDesc.add(generic).add(config);
  opts::variables_map vm;

  try {
    opts::store(opts::command_line_parser(argc, argv).options(cmdlineOptions)
        .positional(positionalOptions)
        .run(), vm);
    opts::notify(vm);
  } catch (exception e) {
    LOG(ERROR) << e.what() << std::endl;
    std::cout << visibleDesc << "\n";
    exit(1);
  }

  if (vm.count("help")) {
    std::cout << "transitmapper [options] <input-feed>\n"
      << VERSION_STR << "\n\n"
      << visibleDesc << "\n";
    exit(0);
  }

  if (vm.count("version")) {
    std::cout << "\n" << VERSION_STR << "\n";
    exit(0);
  }
}
