// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>

#include <exception>
#include <iostream>
#include <string>

#include "octi/_config.h"
#include "octi/config/ConfigReader.h"
#include "util/log/Log.h"

using octi::config::ConfigReader;

using octi::basegraph::BaseGraphType;
using octi::config::OrderMethod;
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
  std::cout << std::setfill(' ') << std::left << "octi (part of LOOM) "
            << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__ << ")"
            << "\n\n(C) 2017-" << YEAR << " " << COPY << "\n"
            << "Authors: " << AUTHORS << "\n\n"
            << "Usage: " << bin << " < graph.json\n\n"
            << "Allowed options:\n\n"
            << "General:\n"
            << std::setw(39) << "  -v [ --version ]"
            << "print version\n"
            << std::setw(39) << "  -h [ --help ]"
            << "show this help message\n"
            << std::setw(39) << "  -m [ --optim-mode ] arg (=heur)"
            << "optimization mode, 'heur' or 'ilp'\n"
            << std::setw(39) << "  --obstacles arg"
            << "GeoJSON file containing obstacle polygons\n"
            << std::setw(39) << "  -g [ --grid-size ] arg (=100%)"
            << "grid cell length, either exact or a\n"
            << std::setw(39) << " "
            << " percentage of input adjacent station distance\n"
            << std::setw(39) << "  -b [ -base-graph ] arg (=octilinear)"
            << "base graph, either ortholinear, octilinear,\n"
            << std::setw(39) << " "
            << " orthoradial, quadtree, octihanan\n\n"
            << "Misc:\n"
            << std::setw(39) << "  --retry-on-error"
            << "retry 85\% of grid size on error, 30 times\n"
            << std::setw(39) << "  --skip-on-error"
            << "skip graph on error\n"
            << std::setw(39) << "  --ilp-num-threads arg (=0)"
            << "number of threads to use by ILP solver,\n"
            << std::setw(39) << " "
            << " 0 means solver default\n"
            << std::setw(39) << "  --hanan-iters arg (=1)"
            << "number of Hanan grid iterations\n"
            << std::setw(39) << "  --loc-search-max-iters arg (=100)"
            << "max local search iterations\n"
            << std::setw(39) << "  --ilp-cache-threshold arg (=inf)"
            << "ILP solve cache treshold\n"
            << std::setw(39) << "  --ilp-time-limit arg (=60)"
            << "ILP solve time limit (seconds), -1 for infinite\n"
            << std::setw(39) << "  --ilp-cache-dir arg (=.)"
            << "ILP cache dir\n"
            << std::setw(39) << "  --ilp-solver arg (=gurobi)"
            << "Preferred ILP solver, either glpk, cbc, or gurobi,\n"
            << std::setw(39) << " "
            << " will fall back if not available.\n"
            << std::setw(39) << "  --write-stats"
            << "write stats to output graph\n"
            << std::setw(39) << "  -D [ --from-dot ]"
            << "input is in dot format\n"
            << std::setw(39) << "  --no-deg2-heur"
            << "don't contract degree 2 nodes\n"
            << std::setw(39) << "  --geo-pen arg (=0)"
            << "enforces lines to follow input geo course\n"
            << std::setw(39) << "  --max-grid-dist arg (=3)"
            << "max grid distance for station candidates\n"
            << std::setw(39) << "  --restr-loc-search"
            << "restrict local search to max grid distance\n"
            << std::setw(39) << "  --edge-order arg (=all)"
            << "method used for initial edge ordering for heur,\n"
            << std::setw(39) << " "
            << " one of num-lines, length, adj-nd-deg,\n"
            << std::setw(39) << " "
            << " adj-nd-ldeg, growth-deg, growth-ldef, all\n"
            << std::setw(39) << "  --density-pen arg (=10)"
            << "restrict local search to max grid distance\n"
            << std::setw(39) << "  --vert-pen arg (=0)"
            << "penalty for vertical edges\n"
            << std::setw(39) << "  --hori-pen arg (=0)"
            << "penalty for horizontal edges\n"
            << std::setw(39) << "  --diag-pen arg (=.5)"
            << "penalty for diagonal edges\n"
            << std::setw(39) << "  --pen-180 arg (=0)"
            << "penalty for 180 deg bends\n"
            << std::setw(39) << "  --pen-135 arg (=1)"
            << "penalty for 135 deg bends\n"
            << std::setw(39) << "  --pen-90 arg (=1.5)"
            << "penalty for 90 deg bends\n"
            << std::setw(39) << "  --pen-45 arg (=2)"
            << "penalty for 45 deg bends\n"
            << std::setw(39) << "  --nd-move-pen arg (=.5)"
            << "penalty for node movement\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) const {
  std::string VERSION_STR = " - unversioned - ";
  std::string baseGraphStr = "octilinear";
  std::string edgeOrderMethod = "all";

  struct option ops[] = {{"version", no_argument, 0, 'v'},
                         {"help", no_argument, 0, 'h'},
                         {"optim-mode", required_argument, 0, 'm'},
                         {"ilp-num-threads", required_argument, 0, 1},
                         {"hanan-iters", required_argument, 0, 2},
                         {"loc-search-max-iters", required_argument, 0, 3},
                         {"ilp-cache-threshold", required_argument, 0, 4},
                         {"ilp-time-limit", required_argument, 0, 5},
                         {"ilp-cache-dir", required_argument, 0, 6},
                         {"ilp-solver", required_argument, 0, 7},
                         {"write-stats", no_argument, 0, 8},
                         {"obstacles", required_argument, 0, 9},
                         {"from-dot", required_argument, 0, 'D'},
                         {"no-deg2-heur", no_argument, 0, 10},
                         {"geo-pen", required_argument, 0, 11},
                         {"max-grid-dist", required_argument, 0, 12},
                         {"restr-loc-search", no_argument, 0, 13},
                         {"edge-order", required_argument, 0, 14},
                         {"density-pen", required_argument, 0, 15},
                         {"grid-size", required_argument, 0, 'g'},
                         {"base-graph", required_argument, 0, 'b'},
                         {"vert-pen", required_argument, 0, 17},
                         {"hori-pen", required_argument, 0, 18},
                         {"diag-pen", required_argument, 0, 19},
                         {"pen-180", required_argument, 0, 20},
                         {"pen-135", required_argument, 0, 21},
                         {"pen-90", required_argument, 0, 22},
                         {"pen-45", required_argument, 0, 23},
                         {"nd-move-pen", required_argument, 0, 24},
                         {"skip-on-error", no_argument, 0, 25},
                         {"retry-on-error", no_argument, 0, 26},
                         {"abort-after", required_argument, 0, 'a'},
                         {0, 0, 0, 0}};

  int c;

  while ((c = getopt_long(argc, argv, ":hvm:Dg:b:", ops, 0)) != -1) {
    switch (c) {
      case 'a':
        cfg->abortAfter = atoi(optarg);
        break;
      case 'D':
        cfg->fromDot = true;
        break;
      case 'm':
        cfg->optMode = optarg;
        break;
      case 1:
        cfg->ilpNumThreads = atoi(optarg);
        break;
      case 2:
        cfg->hananIters = atoi(optarg);
        break;
      case 3:
        cfg->heurLocSearchIters = atoi(optarg);
        break;
      case 4:
        cfg->ilpCacheThreshold = atof(optarg);
        break;
      case 5:
        cfg->ilpTimeLimit = atoi(optarg);
        break;
      case 6:
        cfg->ilpCacheDir = optarg;
        break;
      case 7:
        cfg->ilpSolver = optarg;
        break;
      case 8:
        cfg->writeStats = true;
        break;
      case 9:
        cfg->obstaclePath = optarg;
        break;
      case 10:
        cfg->deg2Heur = false;
        break;
      case 11:
        cfg->enfGeoPen = atof(optarg);
        break;
      case 12:
        cfg->maxGrDist = atof(optarg);
        break;
      case 13:
        cfg->restrLocSearch = true;
        break;
      case 14:
        edgeOrderMethod = optarg;
        break;
      case 15:
        cfg->pens.densityPen = atof(optarg);
        break;
      case 'b':
        baseGraphStr = optarg;
        break;
      case 17:
        cfg->pens.verticalPen = atof(optarg);
        break;
      case 18:
        cfg->pens.horizontalPen = atof(optarg);
        break;
      case 19:
        cfg->pens.diagonalPen = atof(optarg);
        break;
      case 20:
        cfg->pens.p_0 = atof(optarg);
        break;
      case 21:
        cfg->pens.p_135 = atof(optarg);
        break;
      case 22:
        cfg->pens.p_90 = atof(optarg);
        break;
      case 23:
        cfg->pens.p_45 = atof(optarg);
        break;
      case 24:
        cfg->pens.ndMovePen = atof(optarg);
        break;
      case 25:
        cfg->skipOnError = true;
        break;
      case 26:
        cfg->retryOnError = true;
        break;
      case 'g':
        cfg->gridSize = optarg;
        break;
      case 'v':
        std::cout << "octi - (LOOM " << VERSION_FULL << ")" << std::endl;
        exit(0);
      case 'h':
        help(argv[0]);
        exit(0);
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

  if (edgeOrderMethod == "num-lines") {
    cfg->orderMethod = OrderMethod::NUM_LINES;
  } else if (edgeOrderMethod == "length") {
    cfg->orderMethod = OrderMethod::LENGTH;
  } else if (edgeOrderMethod == "adj-nd-deg") {
    cfg->orderMethod = OrderMethod::ADJ_ND_DEGREE;
  } else if (edgeOrderMethod == "adj-nd-ldeg") {
    cfg->orderMethod = OrderMethod::ADJ_ND_LDEGREE;
  } else if (edgeOrderMethod == "growth-deg") {
    cfg->orderMethod = OrderMethod::GROWTH_DEG;
  } else if (edgeOrderMethod == "growth-ldeg") {
    cfg->orderMethod = OrderMethod::GROWTH_LDEG;
  } else if (edgeOrderMethod == "all") {
    cfg->orderMethod = OrderMethod::ALL;
  } else {
    LOG(ERROR) << "Unknown order method " << edgeOrderMethod
               << ", must be one of {num-lines, length, adj-nd-deg, "
                  "adj-nd-ldeg, growth-deg, growth-ldef, all}";
    exit(0);
  }

  if (baseGraphStr == "ortholinear") {
    cfg->baseGraphType = BaseGraphType::GRID;
  } else if (baseGraphStr == "octilinear") {
    cfg->baseGraphType = BaseGraphType::OCTIGRID;
  } else if (baseGraphStr == "hexalinear") {
    cfg->baseGraphType = BaseGraphType::HEXGRID;
  } else if (baseGraphStr == "chulloctilinear") {
    cfg->baseGraphType = BaseGraphType::CONVEXHULLOCTIGRID;
  } else if (baseGraphStr == "porthoradial") {
    cfg->baseGraphType = BaseGraphType::PSEUDOORTHORADIAL;
  } else if (baseGraphStr == "orthoradial") {
    cfg->baseGraphType = BaseGraphType::PSEUDOORTHORADIAL;
  } else if (baseGraphStr == "pseudoorthoradial") {
    cfg->baseGraphType = BaseGraphType::PSEUDOORTHORADIAL;
  } else if (baseGraphStr == "quadtree") {
    cfg->baseGraphType = BaseGraphType::OCTIQUADTREE;
  } else if (baseGraphStr == "octihanan") {
    cfg->baseGraphType = BaseGraphType::OCTIHANANGRID;
  } else if (baseGraphStr == "octihanan") {
    cfg->baseGraphType = BaseGraphType::OCTIHANANGRID;
  } else {
    LOG(ERROR) << "Unknown base graph type " << baseGraphStr << std::endl;
    exit(0);
  }
}
