// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include "lest.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/graph/DirGraph.h"
#include "util/graph/UndirGraph.h"
#include "util/graph/Dijkstra.h"
#include "util/graph/EDijkstra.h"
#include "util/geo/Grid.h"
#include "util/Misc.h"

using lest::approx;
using namespace util;
using namespace util::geo;
using namespace util::graph;

// define LEST cases
const lest::test specification[] = {

// ___________________________________________________________________________
CASE("atof") {
  EXPECT(util::atof("45.534215") == approx(45.534215));
  EXPECT(util::atof("5.534") == approx(5.534));
  EXPECT(util::atof("534") == approx(534));
  EXPECT(util::atof("-534") == approx(-534));
  EXPECT(util::atof("-45.534215") == approx(-45.534215));
  EXPECT(util::atof("-45.534215", 2) == approx(-45.53));


  // TODO: more test cases
}



};

// _____________________________________________________________________________
int main(int argc, char** argv) {
  return(lest::run(specification, argc, argv));
}
