// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include "lest.h"

using lest::approx;
// define LEST cases
const lest::test specification[] = {

// ___________________________________________________________________________
CASE("test") {
  EXPECT(true);
}
};

// _____________________________________________________________________________
int main(int argc, char** argv) {
  return(lest::run(specification, argc, argv));
}
