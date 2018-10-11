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

class HierarchOrderingConfig
    : public std::map<const Edge*, std::map<size_t, Ordering>> {
 public:
  void writeFlatCfg(OrderingConfig* c) const {
    for (auto kv : *this) {
      for (auto ordering : kv.second) {
        (*c)[kv.first].insert((*c)[kv.first].begin(), ordering.second.begin(),
                              ordering.second.end());
      }
    }
  }
};
}
}

#endif  // TRANSITMAP_GRAPH_ORDERINGCONFIG_H_
