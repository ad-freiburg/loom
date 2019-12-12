// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_GRAPH_ROUTE_H_
#define TRANSITMAP_GRAPH_ROUTE_H_

#include <string>
#include <vector>

namespace transitmapper {
namespace graph {

class Route {
 public:
  Route(const std::string& id, const std::string& label, const std::string& color)
  : _id(id), _label(label), _color(color), _orderRelativeTo(0), _numPartners(0) {}

  const std::string& getId() const;
  const std::string& getLabel() const;
  const std::string& getColor() const;

  void setRelativeTo(const Route* r);
  const Route* relativeTo() const;

  size_t getNumCollapsedPartners() const;
  void setNumCollapsedPartners(size_t n);

 private:
  std::string _id, _label, _color;
  const Route* _orderRelativeTo;

  size_t _numPartners;
};

}}

#endif  // TRANSITMAP_GRAPH_ROUTE_H_
