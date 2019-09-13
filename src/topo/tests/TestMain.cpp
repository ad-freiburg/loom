// Copyright 2016
// Author: Patrick Brosi

#include "topo/tests/ContractTest.h"
#include "topo/tests/ContractTest2.h"
#include "topo/tests/TopologicalTest.h"
#include "topo/tests/TopologicalTest2.h"
#include "topo/tests/TopologicalTest3.h"
#include "topo/tests/RestrInfTest.h"

#include "util/Misc.h"

#define TEST(s, o, e) if (!(s o e)) {  std::stringstream ss;  ss << "\n" << __FILE__ << ":" << __LINE__ << ": Test failed!\n  Expected " << #s << " " << #o << " <" << std::to_string(e) << ">, got <" << std::to_string(s) << ">";  throw std::runtime_error(ss.str());}

// _____________________________________________________________________________
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);
  ContractTest2 ct2;
  ContractTest ct;
  TopologicalTest tt;
  TopologicalTest2 tt2;
  TopologicalTest3 tt3;
  RestrInfTest rt;

  TEST(1, ==, util::approx(2.0));

  rt.run();
  ct2.run();
  ct.run();
  tt.run();
  tt2.run();
  tt3.run();
}
