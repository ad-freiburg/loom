// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/graph/Route.h"

using namespace transitmapper;
using namespace graph;

// _____________________________________________________________________________
const std::string& Route::getId() const {
  return _id;
}

// _____________________________________________________________________________
const std::string& Route::getLabel() const {
  return _label;
}

// _____________________________________________________________________________
const std::string& Route::getColor() const {
  return _color;
}

// _____________________________________________________________________________
void Route::setRelativeTo(const Route* r) {
  _orderRelativeTo = r;
}

// _____________________________________________________________________________
const Route* Route::relativeTo() const {
   return _orderRelativeTo;
}
