// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/style/LineStyle.h"

using shared::style::LineStyle;

// _____________________________________________________________________________
void LineStyle::setOutlineCss(const std::string& css) {
  _oCss = css;
}

// _____________________________________________________________________________
const std::string& LineStyle::getOutlineCss() const {
  return _oCss;
}

// _____________________________________________________________________________
void LineStyle::setCss(const std::string& css) {
  _css = css;
}

// _____________________________________________________________________________
const std::string& LineStyle::getCss() const {
  return _css;
}
