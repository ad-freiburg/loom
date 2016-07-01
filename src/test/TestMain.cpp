// Copyright 2016
// Author: Patrick Brosi
//
#include "./lest.h"

using lest::approx;
// using namespace Transitmap;

// define LEST cases
const lest::test specification[] = {

// ___________________________________________________________________________
CASE("dummy") {
  EXPECT(true);
}

};

// _____________________________________________________________________________
int main(int argc, char** argv) {
  return(lest::run(specification, argc, argv));
}
