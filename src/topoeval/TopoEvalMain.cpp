// Copyright 2022
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <stdio.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include "shared/linegraph/LineGraph.h"
#include "util/log/Log.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  shared::linegraph::LineGraph gtGraph;
  shared::linegraph::LineGraph testGraph;

  std::ifstream ifs;

  ifs.open(argv[1]);
  gtGraph.readFromJson(&ifs, false);

  ifs.open(argv[2]);
  testGraph.readFromJson(&ifs, false);

  LOG(INFO) << "Ground truth graph: " << gtGraph.getNds()->size() << " nodes";
  LOG(INFO) << "Test graph: " << testGraph.getNds()->size() << " nodes";

  return (0);
}
