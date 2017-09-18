// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <boost/program_options.hpp>
#include <exception>
#include <iostream>
#include <string>
#include "gtfs2topo/config/ConfigReader.h"
#include "util/log/Log.h"

using gtfs2topo::config::ConfigReader;
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
      "projection",
      opts::value<std::string>(&(cfg->projectionString))
          ->default_value("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 "
                          "+lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m "
                          "+nadgrids=@null +wktext  +no_defs"),
      "map projection as proj4 string")(
      "ignore-gtfs-distances",
      opts::bool_switch(&(cfg->ignoreGtfsDistances))->default_value(false),
      "ignore the distance values in GTFS feed, useful for buggy feeds")(
      "station-aggregation-level",
      opts::value<size_t>(&(cfg->stationAggrLevel))->default_value(2),
      "2 = aggregate based on distance, 1 = aggregate based on feed, 0 = no "
      "aggr")(
      "max-station-aggr-distance",
      opts::value<double>(&(cfg->stationAggrDistance))->default_value(100),
      "maximum aggregation distance between stops (if aggr level > 1)")(
      "max-aggr-distance,d",
      opts::value<double>(&(cfg->maxAggrDistance))->default_value(30),
      "maximum aggregation distance between shared segments")(
      "use-mots,m", opts::value<uint8_t>(&(cfg->useMots))->default_value(255),
      "used mots, as a bitmap (see GTFS reference, 0000001 = 1 means 'only "
      "tram', 0000010 = 2 means 'only subway', 0000100 = 4 means 'only rail', "
      "0001000 = 8 means 'only bus', 0001001 = 3 means 'tram and bus' etc");

  opts::options_description positional("Positional arguments");
  positional.add_options()("input-feed",
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
    opts::store(opts::command_line_parser(argc, argv)
                    .options(cmdlineOptions)
                    .positional(positionalOptions)
                    .run(),
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
