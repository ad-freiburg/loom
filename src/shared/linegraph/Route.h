// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_ROUTE_H_
#define SHARED_LINEGRAPH_ROUTE_H_

#include <string>
#include <vector>

namespace shared {
namespace linegraph {

class Route {
 public:
  Route(const std::string& id, const std::string& label,
        const std::string& color)
      : _id(id),
        _label(label),
        _color(color) {}

  const std::string& getId() const;
  const std::string& getLabel() const;
  const std::string& getColor() const;

 private:
  std::string _id, _label, _color;
};
}
}

#endif  // SHARED_LINEGRAPH_ROUTE_H_
