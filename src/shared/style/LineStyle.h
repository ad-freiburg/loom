// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef SHARED_STYLE_LINESTYLE_H_
#define SHARED_STYLE_LINESTYLE_H_

#include <string>
#include <vector>

namespace shared {
namespace style {

class LineStyle {
 public:
  LineStyle(){};

  const std::vector<double>& getDashArray() const;
  std::string getDashArrayString() const;
  void setDashArray(const std::vector<double>& arr);
  void setDashArray(const std::string& doubleArrayString);

  void setCss(const std::string& css);
  const std::string& getCss() const;

  bool operator==(const LineStyle& other) const {
    return _dashArray == other.getDashArray() && _css == other.getCss();
  }

 private:
  std::vector<double> _dashArray;
  std::string _css;
};
}
}

#endif  // SHARED_STYLE_LINESTYLE_H_
