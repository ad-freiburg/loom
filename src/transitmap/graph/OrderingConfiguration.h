// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_ORDERINGCONFIGURATION_H_
#define TRANSITMAP_GRAPH_ORDERINGCONFIGURATION_H_

#include "./EdgeTripGeom.h"

using std::exception;
using std::string;

namespace transitmapper {
namespace graph {

typedef std::vector<size_t> Ordering;
typedef std::map<EdgeTripGeom*, Ordering> Configuration;

}}

#endif  // TRANSITMAP_GRAPH_ORDERINGCONFIGURATION_H_
