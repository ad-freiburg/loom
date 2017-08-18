// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_ORDERINGCONFIG_H_
#define TRANSITMAP_GRAPH_ORDERINGCONFIG_H_

#include <map>
#include <vector>

namespace transitmapper {
namespace graph {

class Edge;

typedef std::vector<size_t> Ordering;
typedef std::map<const Edge*, Ordering> OrderingConfig;

}}

#endif  // TRANSITMAP_GRAPH_ORDERINGCONFIG_H_
