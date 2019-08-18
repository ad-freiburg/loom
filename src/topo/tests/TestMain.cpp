// Copyright 2016
// Author: Patrick Brosi


#include <cassert>
#include <string>

#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/output/GeoGraphJsonOutput.h"

#define private public
#include "topo/builder/Builder.h"

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
int main(int argc, char** argv) {
  UNUSED(argc);
  UNUSED(argv);

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c <--- a ---> b
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    assert(tg.getNds()->size() == 3);

    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});
    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    ab->pl().addRoute(&l3, 0);

    ac->pl().addRoute(&l1, 0);
    ac->pl().addRoute(&l2, 0);
    ac->pl().addRoute(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::Builder builder(&cfg);

    builder.combineEdges(ab, ac, a, &tg);

    assert(tg.getNds()->size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    assert(b->getAdjList().front() == c->getAdjList().front());
    assert(b->getAdjList().front()->pl().getRoutes().size() == 3);

    // check geometry orientation
    assert(b->getAdjList().front()->pl().getGeom()->front().getX() == approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(b->getAdjList().front()->pl().getGeom()->back().getX() == approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getRoutes()) {
      assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c <--- a <--- b
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    assert(tg.getNds()->size() == 3);

    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ba->pl().addRoute(&l1, 0);
    ba->pl().addRoute(&l2, 0);
    ba->pl().addRoute(&l3, 0);

    ac->pl().addRoute(&l1, 0);
    ac->pl().addRoute(&l2, 0);
    ac->pl().addRoute(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::Builder builder(&cfg);

    builder.combineEdges(ba, ac, a, &tg);

    util::geo::output::GeoGraphJsonOutput out;
    out.print(tg, std::cout);

    assert(tg.getNds()->size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    assert(b->getAdjList().front() == c->getAdjList().front());
    assert(b->getAdjList().front()->pl().getRoutes().size() == 3);

    // check geometry orientation
    assert(b->getAdjList().front()->pl().getGeom()->front().getX() == approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(b->getAdjList().front()->pl().getGeom()->back().getX() == approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getRoutes()) {
      assert(r.direction == 0);
    }

  }
}
