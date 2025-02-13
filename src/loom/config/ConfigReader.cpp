
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>
#include <exception>
#include <iostream>
#include <string>
#include "loom/_config.h"
#include "loom/config/ConfigReader.h"
#include "loom/config/LoomConfig.h"
#include "util/log/Log.h"

using loom::config::ConfigReader;

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
  std::cout << std::setfill(' ') << std::left << "loom (part of LOOM) "
            << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__ << ")"
            << "\n\n(C) 2017-" << YEAR << " " << COPY << "\n"
            << "Authors: " << AUTHORS << "\n\n"
            << "Usage: " << bin << " < graph.json\n\n"
            << "Allowed options:\n\n"
            << "General:\n"
            << std::setw(43) << "  -v [ --version ]"
            << "print version\n"
            << std::setw(43) << "  -h [ --help ]"
            << "show this help message\n"
            << std::setw(43) << "  --no-untangle"
            << "Don't apply untangling rules\n"
            << std::setw(43) << "  --no-prune"
            << "Don't apply pruning rules\n"
            << std::setw(43) << "  -m [ --optim-method ] arg (=comb-no-ilp)"
            << "Optimization method, one of ilp-naive, ilp,\n"
            << std::setw(43) << " "
            << " comb, comb-no-ilp, exhaust, hillc, hillc-random, anneal,\n"
            << std::setw(43) << " "
            << " anneal-random, greedy, greedy-lookahead, null\n"
            << std::setw(43) << "  --same-seg-cross-pen arg (=4)"
            << "Penalty for same-segment crossings\n"
            << std::setw(43) << "  --diff-seg-cross-pen arg (=1)"
            << "Penalty for diff-segment crossings\n"
            << std::setw(43) << "  --in-stat-cross-pen-same-seg arg (=12)"
            << "Penalty for same-segment crossings at stations\n"
            << std::setw(43) << "  --in-stat-cross-pen-diff-seg arg (=3)"
            << "Penalty for diff-segment crossings at stations\n"
            << std::setw(43) << "  --sep-pen arg (=3)"
            << "Penalty for separations\n"
            << std::setw(43) << "  --in-stat-sep-pen arg (=9)"
            << "Penalty for separations at stations\n\n"
            << "Misc:\n"
            << std::setw(43) << "  -D [ --from-dot ]"
            << "input is in dot format\n"
            << std::setw(43) << "  --output-stats"
            << "Print stats to stdout\n"
            << std::setw(43) << "  --write-stats"
            << "Write stats to output\n"
            << std::setw(43) << "  --ilp-solver arg (=gurobi)"
            << "Preferred ILP solver, either glpk, cbc, or gurobi.\n"
            << std::setw(43) << " "
            << "Will fall back if not available.\n"
            << std::setw(43) << "  --ilp-num-threads arg (=0)"
            << "Number of threads to use by ILP solver,\n"
            << std::setw(43) << " "
            << " 0 means solver default\n"
            << std::setw(43) << "  --ilp-time-limit arg (=-1)"
            << "ILP solve time limit (seconds), -1 for infinite\n"
            << std::setw(43) << "  --dbg-output-path arg (=.)"
            << "Path used for debug output\n"
            << std::setw(43) << "  --output-optgraph"
            << "Output optimization graph to debug path\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) const {
  struct option ops[] = {
      {"version", no_argument, 0, 'v'},
      {"help", no_argument, 0, 'h'},
      {"no-untangle", no_argument, 0, 1},
      {"output-stats", no_argument, 0, 2},
      {"no-prune", no_argument, 0, 3},
      {"same-seg-cross-pen", required_argument, 0, 4},
      {"diff-seg-cross-pen", required_argument, 0, 9},
      {"sep-pen", required_argument, 0, 5},
      {"from-dot", required_argument, 0, 'D'},
      {"in-stat-sep-pen", required_argument, 0, 8},
      {"in-stat-cross-pen-same-seg", required_argument, 0, 6},
      {"in-stat-cross-pen-diff-seg", required_argument, 0, 7},
      {"ilp-num-threads", required_argument, 0, 10},
      {"ilp-time-limit", required_argument, 0, 11},
      {"ilp-solver", required_argument, 0, 12},
      {"optim-method", required_argument, 0, 'm'},
      {"optim-runs", required_argument, 0, 13},
      {"dbg-output-path", required_argument, 0, 14},
      {"output-optgraph", required_argument, 0, 15},
      {"write-stats", no_argument, 0, 16},
      {0, 0, 0, 0}};

  int c;
  while ((c = getopt_long(argc, argv, ":hvm:D", ops, 0)) != -1) {
    switch (c) {
      case 'h':
        help(argv[0]);
        exit(0);
      case 'v':
        std::cout << "loom - (LOOM " << VERSION_FULL << ")" << std::endl;
        exit(0);
      case 1:
        cfg->untangleGraph = false;
        break;
      case 2:
        cfg->outputStats = true;
        break;
      case 3:
        cfg->pruneGraph = false;
        break;
      case 'm':
        cfg->optimMethod = optarg;
        break;
      case 4:
        cfg->crossPenMultiSameSeg = atof(optarg);
        break;
      case 5:
        cfg->separationPenWeight = atof(optarg);
        break;
      case 6:
        cfg->stationCrossWeightSameSeg = atof(optarg);
        break;
      case 7:
        cfg->stationCrossWeightDiffSeg = atof(optarg);
        break;
      case 8:
        cfg->stationSeparationWeight = atof(optarg);
        break;
      case 9:
        cfg->crossPenMultiDiffSeg = atof(optarg);
        break;
      case 10:
        cfg->ilpNumThreads = atoi(optarg);
        break;
      case 11:
        cfg->ilpTimeLimit = atof(optarg);
        break;
      case 12:
        cfg->ilpSolver = optarg;
        break;
      case 13:
        cfg->optimRuns = atoi(optarg);
        break;
      case 14:
        cfg->dbgPath = optarg;
        break;
      case 15:
        cfg->outOptGraph = true;
        break;
      case 16:
        cfg->writeStats = true;
        break;
      case 'D':
        cfg->fromDot = true;
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
}
