// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>
//
#include "util/geo/output/GeoJsonOutput.h"

using namespace util;
using namespace geo;
using namespace output;

// _____________________________________________________________________________
GeoJsonOutput::GeoJsonOutput(std::ostream& str) : GeoJsonOutput(str, false) {}

// _____________________________________________________________________________
GeoJsonOutput::GeoJsonOutput(std::ostream& str, bool raw) : _wr(&str, 10, true) {
  if (!raw) {
    _wr.obj();
    _wr.keyVal("type", "FeatureCollection");
    _wr.key("features");
    _wr.arr();
  }
}

// _____________________________________________________________________________
GeoJsonOutput::GeoJsonOutput(std::ostream& str, json::Val attrs)
    : _wr(&str, 10, true) {
  _wr.obj();
  _wr.keyVal("type", "FeatureCollection");
  _wr.key("properties");
  _wr.val(attrs);
  _wr.key("features");
  _wr.arr();
}

// _____________________________________________________________________________
GeoJsonOutput::~GeoJsonOutput() { flush(); }

// _____________________________________________________________________________
void GeoJsonOutput::flush() { _wr.closeAll(); }
