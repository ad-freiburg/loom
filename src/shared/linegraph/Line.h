// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_LINEGRAPH_LINE_H_
#define SHARED_LINEGRAPH_LINE_H_

#include <string>
#include <vector>

namespace shared {
namespace linegraph {

class Line {
 public:
  Line(const std::string& id, const std::string& label,
       const std::string& color)
      : _id(id), _label(label), _color(color) {}

  const std::string& id() const;
  const std::string& label() const;
  const std::string& color() const;
  void setColor(const std::string& c) { _color = c; };

 private:
  std::string _id, _label, _color;
};
}
}

#endif  // SHARED_LINEGRAPH_LINE_H_
