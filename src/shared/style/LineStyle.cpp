// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <sstream>
#include "shared/style/LineStyle.h"

using shared::style::LineStyle;

// ____________________________________________________________________________
const std::vector<double>& LineStyle::getDashArray() const {
  return _dashArray;
}

// ____________________________________________________________________________
std::string LineStyle::getDashArrayString() const {
  std::stringstream ss;
  for (auto d : _dashArray) ss << d << " ";

  return ss.str();
}

// ____________________________________________________________________________
void LineStyle::setCss(const std::string& css) { _css = css; }

// ____________________________________________________________________________
const std::string& LineStyle::getCss() const { return _css; }

// ____________________________________________________________________________
void LineStyle::setDashArray(const std::vector<double>& arr) {
  _dashArray = arr;
}

// ____________________________________________________________________________
void LineStyle::setDashArray(const std::string& doubleArrayString) {
  _dashArray.clear();
  std::stringstream ss(doubleArrayString);
  for (double a; ss >> a;) _dashArray.push_back(a);
}
