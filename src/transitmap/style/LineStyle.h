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
  void setDashArray(const std::vector<double>& arr);
  void setDashArray(const std::string& doubleArrayString);

 private:
  std::vector<double> _dashArray;

};

}}


#endif  // TRANSITMAP_STYLE_LINESTYLE_H_