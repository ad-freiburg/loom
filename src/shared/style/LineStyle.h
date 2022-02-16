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

  void setOutlineCss(const std::string& css);
  const std::string& getOutlineCss() const;

  void setCss(const std::string& css);
  const std::string& getCss() const;

 private:
  std::string _css, _oCss;
};
}
}

#endif  // SHARED_STYLE_LINESTYLE_H_
