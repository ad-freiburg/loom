// Copyright 2016
// Author: Patrick Brosi


#include <cassert>
#include <string>

#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "topo/tests/ContractTest3.h"
#include "topo/tests/TopoTestUtil.h"
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
void ContractTest3::run() {

  // ___________________________________________________________________________
  {
    /*
     * node contraction of c, b
     *      1      1      1
     *  a ---> b ---> c ----> e
     *         ^ -> /
     *         | 1 / 1 ->
     *         d -/
     *
     */

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});
    auto e = tg.addNd({{15.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});
    auto ce = tg.addEdg(c, e, {{{10.0, 0.0}, {15.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, d);
    dc->pl().addRoute(&l1, d);
    ce->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(c, d, &tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    // assert(tg.getEdg(a, b)->pl().getRoutes().size() == 1);
    // assert(tg.getEdg(d, b)->pl().getRoutes().size() == 1);

    // assert(b->pl().connOccurs(&l1, tg.getEdg(a, b), tg.getEdg(b, d)));
    // assert(validExceptions(&tg));

    // exit(1);
  }
  // ___________________________________________________________________________
  {
    /*
     *    1      1
     * a ---> b ---> c
     *        ^
     *        | 1
     *        |
     *        d
     */

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto db = tg.addEdg(d, b, {{{5.0, -5.0}, {5.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    db->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    b->pl().addConnExc(&l1, bc, db);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

    assert(!d->pl().connOccurs(&l1, tg.getEdg(a, d), tg.getEdg(d, c)));
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    /*
     *        e
     *        ^
     *        | 1
     *    1   |  1
     * a ---> b ---> c
     *        ^
     *        | 1
     *        |
     *        d
     */

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});
    auto e = tg.addNd({{5.0, 5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{5.0, 0.0}, {5.0, 5.0}}});
    auto db = tg.addEdg(d, b, {{{5.0, -5.0}, {5.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    db->pl().addRoute(&l1, 0);
    be->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    b->pl().addConnExc(&l1, bc, db);
    b->pl().addConnExc(&l1, be, db);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

    assert(!d->pl().connOccurs(&l1, tg.getEdg(a, d), tg.getEdg(d, c)));
    // because this connection was already possible before!
    assert(d->pl().connOccurs(&l1, tg.getEdg(e, d), tg.getEdg(d, c)));
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    /*
     *        e
     *        ^
     *        | 1
     *    1   |  1
     * a ---> b ---> c
     *        ^
     *        | 1
     *        |
     *        d
     */

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});
    auto e = tg.addNd({{5.0, 5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{5.0, 0.0}, {5.0, 5.0}}});
    auto db = tg.addEdg(d, b, {{{5.0, -5.0}, {5.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    db->pl().addRoute(&l1, 0);
    be->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    b->pl().addConnExc(&l1, ab, db);
    b->pl().addConnExc(&l1, be, bc);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

    assert(!d->pl().connOccurs(&l1, tg.getEdg(a, d), tg.getEdg(d, c)));
    // because this connection is made possible by bd
    assert(d->pl().connOccurs(&l1, tg.getEdg(e, d), tg.getEdg(d, c)));
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    /*
     *        e
     *        ^
     *        | 1
     *    1   |  1
     * a ---> b ---> c
     *        ^
     *        | 1
     *        |
     *        d
     */

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});
    auto e = tg.addNd({{5.0, 5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{5.0, 0.0}, {5.0, 5.0}}});
    auto db = tg.addEdg(d, b, {{{5.0, -5.0}, {5.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    db->pl().addRoute(&l1, 0);
    be->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    b->pl().addConnExc(&l1, bc, db);
    b->pl().addConnExc(&l1, be, bc);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(!d->pl().connOccurs(&l1, tg.getEdg(a, d), tg.getEdg(d, c)));
    assert(!d->pl().connOccurs(&l1, tg.getEdg(e, d), tg.getEdg(d, c)));
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    /*
     *   1,2     1      1,2
     * a ---> b <--- c <---- e
     *        |      ^
     *        | 2    | 2
     *        v      |
     *        d -----|
     */

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});
    auto e = tg.addNd({{15.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto cb = tg.addEdg(c, b, {{{10.0, 0.0}, {5.0, 0.0}}});
    auto ec = tg.addEdg(e, c, {{{15.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    cb->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l2, 0);
    dc->pl().addRoute(&l2, 0);
    ec->pl().addRoute(&l1, 0);
    ec->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(b, c, &tg);
    builder.combineNodes(c, d, &tg);

    assert(d->pl().connOccurs(&l1, tg.getEdg(a, b), tg.getEdg(d, e)));
    assert(d->pl().connOccurs(&l2, tg.getEdg(a, b), tg.getEdg(d, e)));

    assert(validExceptions(&tg));
  }
}
