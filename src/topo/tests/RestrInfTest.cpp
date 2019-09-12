// Copyright 2016
// Author: Patrick Brosi

#include <cassert>
#include <map>
#include <set>
#include <string>

#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/tests/RestrInfTest.h"
#include "topo/tests/TopoTestUtil.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/output/GeoGraphJsonOutput.h"

#define private public
#include "topo/restr/RestrInferrer.h"

class approx {
 public:
  explicit approx(double magnitude)
      : _epsilon{std::numeric_limits<float>::epsilon() * 100},
        _magnitude{magnitude} {}

  friend bool operator==(double lhs, approx const& rhs) {
    return std::abs(lhs - rhs._magnitude) < rhs._epsilon;
  }

  friend bool operator==(approx const& lhs, double rhs) {
    return operator==(rhs, lhs);
  }
  friend bool operator!=(double lhs, approx const& rhs) {
    return !operator==(lhs, rhs);
  }
  friend bool operator!=(approx const& lhs, double rhs) {
    return !operator==(rhs, lhs);
  }

  friend bool operator<=(double lhs, approx const& rhs) {
    return lhs < rhs._magnitude || lhs == rhs;
  }
  friend bool operator<=(approx const& lhs, double rhs) {
    return lhs._magnitude < rhs || lhs == rhs;
  }
  friend bool operator>=(double lhs, approx const& rhs) {
    return lhs > rhs._magnitude || lhs == rhs;
  }
  friend bool operator>=(approx const& lhs, double rhs) {
    return lhs._magnitude > rhs || lhs == rhs;
  }

 private:
  double _epsilon;
  double _magnitude;
};

// _____________________________________________________________________________
void RestrInfTest::run() {
  // ___________________________________________________________________________
  {
    topo::restr::RestrGraph rg;
    auto a = rg.addNd();
    auto b = rg.addNd();
    auto c = rg.addNd();
    auto d = rg.addNd();

    auto ab = rg.addEdg(a, b, util::geo::PolyLine<double>({0, 10}, {0, 20}));
    auto bc = rg.addEdg(b, c, util::geo::PolyLine<double>({0, 20}, {0, 30}));
    auto cd = rg.addEdg(c, d, util::geo::PolyLine<double>({0, 30}, {0, 40}));

    transitmapper::graph::Route l1("1", "1", "red");

    ab->pl().routes.insert(&l1);
    bc->pl().routes.insert(&l1);
    cd->pl().routes.insert(&l1);

    UNUSED(bc);

    std::set<RestrEdge*> from = {ab};
    std::set<RestrEdge*> to = {cd};

    topo::restr::CostFunc cFunc(&l1, 100);
    double cost = EDijkstra::shortestPath(from, to, cFunc);

    assert(cost == approx(20));
  }
}
