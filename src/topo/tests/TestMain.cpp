// Copyright 2016
// Author: Patrick Brosi

#include "topo/tests/ContractTest.h"
#include "topo/tests/ContractTest2.h"
#include "topo/tests/TopologicalTest.h"
#include "topo/tests/TopologicalTest2.h"
#include "topo/tests/TopologicalTest3.h"
#include "topo/tests/RestrInfTest.h"

// _____________________________________________________________________________
int main(int argc, char** argv) {
  ContractTest2 ct2;
  ContractTest ct;
  TopologicalTest tt;
  TopologicalTest2 tt2;
  TopologicalTest3 tt3;
  RestrInfTest rt;

  rt.run();
  ct2.run();
  ct.run();
  tt.run();
  tt2.run();
  tt3.run();
}
