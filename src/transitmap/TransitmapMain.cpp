// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi 

#include <unistd.h>
#include <iostream>
#include <string>
#include <set>
#include <stdio.h>

// using namespace Transitmap;
using std::string;

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // Disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());

  return(0);
}
