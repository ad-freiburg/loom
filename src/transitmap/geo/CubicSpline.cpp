// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "../util/Geo.h"
#include "./CubicSpline.h"

using namespace transitmapper;
using namespace geo;

// _____________________________________________________________________________
CubicSpline::CubicSpline(const PolyLine& knots) : _didRender(false) {
  findPolynoms(knots.getLine(), &_set);
}

// _____________________________________________________________________________
CubicSpline::CubicSpline(const util::geo::Line& knots) : _didRender(false) {
  findPolynoms(knots, &_set);
}

// _____________________________________________________________________________
double CubicPolynom::valueAt(double atx) const {
  double dx = atx - x;
  return a + b * dx + c * dx * dx + d * dx * dx * dx;
}

// _____________________________________________________________________________
PolyLine CubicSpline::render(double d) const {
  if (!_didRender) {
    // TODO:render

    while (true) {
 

    }

    _didRender = true;
  }
  return _rendered;
}


// _____________________________________________________________________________
void CubicSpline::findPolynoms(const util::geo::Line& knots, std::vector<CubicPolynom>* ret)
const {
  /**
   * based on
   * https://en.wikipedia.org/w/index.php?title=Spline_%28mathematics%29&oldid=288288033#Algorithm_for_computing_natural_cubic_splines
   */

  std::vector<double> a(knots.size());
  for (auto& p : knots) {
    a.push_back(p.get<1>());
  }

  std::vector<double> b(knots.size()),
    d(knots.size()),
    h(knots.size());

  for (size_t i = 0; i < knots.size(); ++i) {
    h.push_back(knots[i+1].get<0>() - knots[i].get<0>());
  }

  std::vector<double> alpha;
  for (size_t i = 0; i < knots.size(); ++i) {
    alpha.push_back(3 * (a[i+1]-a[i])/h[i] - 3 * (a[i]-a[i-1])/h[i-1]);
  }

  std::vector<double>
    c(knots.size()+1),
    l(knots.size()+1),
    mu(knots.size()+1),
    z(knots.size()+1);

  l[0] = 1;
  mu[0] = 0;
  z[0] = 0;

  for (int64_t i = knots.size()-1; i >= 0; --i) {
    c[i] = z[i] - mu[i] * c[i+1];
    b[i] = (a[i+1]-a[i])/h[i]-h[i]*(c[i+1]+2*c[i])/3;
    d[i] = (c[i+1] - c[i]) / 3/h[i];
  }

  for (size_t i = 0; i < knots.size(); ++i) {
    ret->push_back(CubicPolynom(a[i], b[i], c[i], d[i], knots[i].get<0>()));
  }
}
