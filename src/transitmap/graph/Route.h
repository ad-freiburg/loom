// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_ROUTE_H_
#define TRANSITMAP_GRAPH_ROUTE_H_

namespace transitmapper {
namespace graph {

struct Route {
  Route(const std::string& id, std::string& label, const std::string& color)
  : id(id), label(label), color(color) {}

  std::string id, label, color;
};

}}

#endif  // TRANSITMAP_GRAPH_ROUTE_H_