// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "LineStyle.h"

using namespace transitmapper;
using namespace style;

// ____________________________________________________________________________
const std::vector<double>& LineStyle::getDashArray() const {
  return _dashArray;
}

// ____________________________________________________________________________
void LineStyle::setDashArray(const std::vector<double>& arr) {

}

// ____________________________________________________________________________
void LineStyle::setDashArray(const std::string& doubleArrayString) {

}
