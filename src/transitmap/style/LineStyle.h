// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_STYLE_LINESTYLE_H_
#define TRANSITMAP_STYLE_LINESTYLE_H_

#include <vector>
#include <string>

namespace transitmapper {
namespace style {

class LineStyle {

 public:
  LineStyle() {};

  const std::vector<double>& getDashArray() const;
  std::string getDashArrayString() const;
  void setDashArray(const std::vector<double>& arr);
  void setDashArray(const std::string& doubleArrayString);

  void setCss(const std::string& css);
  const std::string& getCss();

 private:
  std::vector<double> _dashArray;
  std::string _css;

};

}}


#endif  // TRANSITMAP_STYLE_LINESTYLE_H_
