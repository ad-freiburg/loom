// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>

#include <exception>
#include <iostream>
#include <string>

#include "ad/cppgtfs/gtfs/flat/Route.h"
#include "gtfs2graph/_config.h"
#include "gtfs2graph/config/ConfigReader.h"
#include "util/String.h"
#include "util/log/Log.h"

using gtfs2graph::config::ConfigReader;

using std::exception;
using std::string;
using std::vector;

static const char* YEAR = &__DATE__[7];
static const char* COPY =
    "University of Freiburg - Chair of Algorithms and Data Structures";
static const char* AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

// _____________________________________________________________________________
ConfigReader::ConfigReader() {}

// _____________________________________________________________________________
void ConfigReader::help(const char* bin) const {
  std::cout
      << std::setfill(' ') << std::left << "gtfs2graph (part of LOOM) "
      << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__ << ")"
      << "\n\n(C) 2017-" << YEAR << " " << COPY << "\n"
      << "Authors: " << AUTHORS << "\n\n"
      << "Usage: " << bin << " <GTFS FEED>\n\n"
      << "Allowed options:\n\n"
      << "General:\n"
      << std::setw(36) << "  -v [ --version ]"
      << "print version\n"
      << std::setw(36) << "  -h [ --help ]"
      << "show this help message\n"
      << std::setw(36) << "  -m [ --mots ] arg (=all)"
      << "MOTs to calculate shapes for, comma sep.,\n"
      << std::setw(36) << " "
      << "  either as string "
         "{all, tram | streetcar,\n"
      << std::setw(36) << " "
      << "  subway | metro, rail | train, bus,\n"
      << std::setw(36) << " "
      << "  ferry | boat | ship, cablecar, gondola,\n"
      << std::setw(36) << " "
      << "  funicular, coach, mono-rail | monorail,\n"
      << std::setw(36) << " "
      << "  trolley | trolleybus | trolley-bus} or\n"
      << std::setw(36) << " "
      << "  as GTFS mot codes\n"
      << std::setw(36) << " "
      << "  funicular, coach} or as GTFS mot codes\n"
      << std::setw(36) << "  -r [ --route ]"
      << "Routes to calculate shapes for, comma separator.\n"
      << std::setw(36) << "  -p [ --prune-threshold ] arg (=0)"
      << "Threshold for pruning of seldomly occuring\n"
      << std::setw(36) << " " << "  lines, between 0 and 1\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) const {
  std::string motStr = "all";
  std::string routeIds = "-1";
  double pruneThreshold = 0;

  struct option ops[] = {{"version", no_argument, 0, 'v'},
                         {"help", no_argument, 0, 'h'},
                         {"mots", required_argument, 0, 'm'},
                         {"route", required_argument, 0, 'r'},
                         {"prune-threshold", required_argument, 0, 'p'},
                         {0, 0, 0, 0}};

  int c;
  while ((c = getopt_long(argc, argv, ":hvim:r:p:", ops, 0)) != -1) {
    switch (c) {
      case 'h':
        help(argv[0]);
        exit(0);
      case 'v':
        std::cout << "gtfs2graph - (LOOM " << VERSION_FULL << ")" << std::endl;
        exit(0);
      case 'm':
        motStr = optarg;
        break;
      case 'r':
        routeIds = optarg;
        break;
      case 'p':
        pruneThreshold = atof(optarg);
        break;
      case ':':
        std::cerr << argv[optind - 1];
        std::cerr << " requires an argument" << std::endl;
        exit(1);
      case '?':
        std::cerr << argv[optind - 1];
        std::cerr << " option unknown" << std::endl;
        exit(1);
        break;
      default:
        std::cerr << "Error while parsing arguments" << std::endl;
        exit(1);
        break;
    }
  }

  if (optind == argc) {
    std::cerr << "No input GTFS feed specified." << std::endl;
    exit(1);
  }

  cfg->inputFeedPath = argv[optind];
  cfg->pruneThreshold = pruneThreshold;

  for (auto sMotStr : util::split(motStr, ',')) {
    for (auto mot :
         ad::cppgtfs::gtfs::flat::Route::getTypesFromString(sMotStr)) {
      cfg->useMots.insert(mot);
    }
  }

  for (auto sRouteStr : util::split(routeIds, ',')) {
    if (sRouteStr == "-1") continue;
    cfg->useRoutes.insert(sRouteStr);
  }
}