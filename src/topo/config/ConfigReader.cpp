// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>

#include <exception>
#include <iostream>
#include <string>

#include "topo/_config.h"
#include "topo/config/ConfigReader.h"
#include "util/log/Log.h"

using topo::config::ConfigReader;

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
  std::cout << std::setfill(' ') << std::left << "topo (part of LOOM) "
            << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__ << ")"
            << "\n\n(C) 2017-" << YEAR << " " << COPY << "\n"
            << "Authors: " << AUTHORS << "\n\n"
            << "Usage: " << bin << " < linegraph.json\n\n"
            << "Allowed options:\n\n"
            << "General:\n"
            << std::setw(40) << "  -v [ --version ]"
            << "print version\n"
            << std::setw(40) << "  -h [ --help ]"
            << "show this help message\n"
            << std::setw(40) << "  -d [ --max-aggr-dist ] arg (=50)"
            << "maximum distance between segments\n"
            << std::setw(40) << "  --write-stats"
            << "write statistics to output file\n"
            << std::setw(40) << "  --no-infer-restrs"
            << "don't infer turn restrictions\n"
            << std::setw(40) << "  --infer-restr-max-dist arg (=[-d])"
            << "max dist for considered edges for turn restrictions\n"
            << std::setw(40) << "  --max-comp-dist arg (=10000)"
            << "max distance between nodes in component, in meters\n"
            << std::setw(40) << "  --sample-dist arg (=5)"
            << "sample length for map construction, in pseudometers\n"
            << std::setw(40) << "  --max-length-dev arg (=500)"
            << "maximum distance deviation for turn restrictions infer\n"
            << std::setw(40) << "  --turn-restr-full-turn-angle arg (=0)"
            << "turn angles smaller than this will count as full turn\n"
            << std::setw(40) << "  --turn-restr-full-turn-pen arg (=0)"
            << "penalty for full turns during turn restriction infer\n"
            << std::setw(40) << "  --random-colors"
            << "fill missing colors with random colors\n"
            << std::setw(40) << "  --write-components"
            << "write graph component ID to edge attributes\n"
            << std::setw(40) << "  --write-components-path"
            << "write graph components as separated files to given path\n"
            << std::setw(40) << "  --smooth (=0)"
            << "smooth output graph edge geometries\n"
            << std::setw(40) << "  --aggr-stats"
            << "aggregate stats with existing from input\n";
}

// _____________________________________________________________________________
void ConfigReader::read(TopoConfig* cfg, int argc, char** argv) const {
  std::string motStr = "all";

  struct option ops[] = {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"max-aggr-dist", required_argument, 0, 'd'},
      {"no-infer-restrs", no_argument, 0, 1},
      {"write-stats", no_argument, 0, 2},
      {"max-length-dev", required_argument, 0, 3},
      {"infer-restr-max-dist", required_argument, 0, 4},
      {"write-components", no_argument, 0, 5},
      {"write-components-path", required_argument, 0, 6},
      {"turn-restr-full-turn-pen", required_argument, 0, 7},
      {"random-colors", no_argument, 0, 8},
      {"sample-dist", required_argument, 0, 9},
      {"max-comp-dist", required_argument, 0, 10},
      {"smooth", required_argument, 0, 11},
      {"turn-restr-full-turn-angle", required_argument, 0, 12},
      {"aggr-stats", no_argument, 0, 13},
      {0, 0, 0, 0}};

  double turnRestrDiff = -1;

  int c;
  while ((c = getopt_long(argc, argv, ":hvd:", ops, 0)) != -1) {
    switch (c) {
      case 'h':
        help(argv[0]);
        exit(0);
      case 'v':
        std::cout << "topo - (LOOM " << VERSION_FULL << ")" << std::endl;
        exit(0);
      case 'd':
        cfg->maxAggrDistance = atof(optarg);
        break;
      case 1:
        cfg->noInferRestrs = true;
        break;
      case 2:
        cfg->outputStats = true;
        break;
      case 3:
        cfg->maxAggrDistance = atof(optarg);
        break;
      case 4:
        turnRestrDiff = atof(optarg);
        break;
      case 5:
        cfg->writeComponents = true;
        break;
      case 6:
        cfg->componentsPath = optarg;
        break;
      case 7:
        cfg->turnInferFullTurnPen = atof(optarg);
        break;
      case 8:
        cfg->randomColors = true;
        break;
      case 9:
        cfg->segmentLength = atof(optarg);
        break;
      case 10:
        cfg->connectedCompDist = atof(optarg);
        break;
      case 11:
        cfg->smooth = atof(optarg);
        break;
      case 12:
        cfg->fullTurnAngle = atof(optarg);
        break;
      case 13:
        cfg->aggregateStats = true;
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

  if (turnRestrDiff >= 0)
    cfg->maxTurnRestrCheckDist = turnRestrDiff;
  else
    cfg->maxTurnRestrCheckDist = cfg->maxAggrDistance;
}
