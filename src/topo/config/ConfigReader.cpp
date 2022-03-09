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
            << "\n\n(C) " << YEAR << " " << COPY << "\n"
            << "Authors: " << AUTHORS << "\n\n"
            << "Usage: " << bin << " < linegraph.json\n\n"
            << "Allowed options:\n\n"
            << "General:\n"
            << std::setw(35) << "  -v [ --version ]"
            << "print version\n"
            << std::setw(35) << "  -h [ --help ]"
            << "show this help message\n"
            << std::setw(35) << "  -d [ --max-aggr-dist ] arg (=40)"
            << "maximum distance between segments\n"
            << std::setw(35) << "  --write-stats"
            << "write statistics to output file\n"
            << std::setw(35) << "  --no-infer-restrs"
            << "don't infer turn restrictions\n"
            << std::setw(35) << "  --max-length-dev arg (=500)"
            << "maxumum distance deviation for turn restrictions infer\n";
}

// _____________________________________________________________________________
void ConfigReader::read(TopoConfig* cfg, int argc, char** argv) const {
  std::string motStr = "all";

  struct option ops[] = {{"version", no_argument, 0, 'v'},
                         {"help", no_argument, 0, 'h'},
                         {"max-aggr-dist", required_argument, 0, 'd'},
                         {"no-infer-restrs", no_argument, 0, 1},
                         {"write-stats", no_argument, 0, 2},
                         {"max-length-dev", required_argument, 0, 3},
                         {0, 0, 0, 0}};

  char c;
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
