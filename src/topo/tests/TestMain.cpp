// Copyright 2016
// Author: Patrick Brosi

#include "topo/tests/ContractTest.h"
#include "topo/tests/TopologicalTest.h"
#include "topo/tests/TopologicalTest2.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  ContractTest ct;
  TopologicalTest tt;
  TopologicalTest2 tt2;
  ct.run();
  tt.run();
  tt2.run();
}
