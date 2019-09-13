// Copyright 2016
// Author: Patrick Brosi


#include <cassert>
#include <string>

#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/tests/TopologicalTest3.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "topo/tests/TopoTestUtil.h"

#define private public
#include "topo/mapconstructor/MapConstructor.h"

using util::approx;

void TopologicalTest3::run() {
  // ___________________________________________________________________________
  {
    //    1, 2           1                    1              1, 2
    // X -----> A --------------> B --------------------> C -------> Y
    //          | ----------------------------------------> D
    //                                2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{250.0, 0.0}});
    auto c = tg.addNd({{500.0, 0}});
    auto d = tg.addNd({{510.0, 5}});
    auto y = tg.addNd({{800.0, 1.0}});
    auto x = tg.addNd({{-300.0, 0.0}});

    auto xa = tg.addEdg(x, a, {{{-300.0, 0.0},  {0.0, 0.0}}});
    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {250.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{250.0, 0.0}, {500.0, 0}}});
    auto ad = tg.addEdg(a, d, {{{0.0, 0.0}, {510.0, 5.0}}});
    auto cy = tg.addEdg(c, y, {{{500.0, 0.0}, {800.0, 1.0}}});


    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    xa->pl().addRoute(&l1, 0);
    xa->pl().addRoute(&l2, 0);

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    ad->pl().addRoute(&l2, 0);
    cy->pl().addRoute(&l1, 0);
    cy->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;
    cfg.minSegLength = 20;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(tg.getNds()->size() == 5);

    a = x->getAdjList().front()->getOtherNd(x);
    d = y->getAdjList().front()->getOtherNd(y);
    for (auto nd : *tg.getNds()) if (nd != a && nd != d && nd != x && nd != y) b = nd;
  }

  // ___________________________________________________________________________
  {
    //    1, 2           1                    1              1, 2
    // X -----> A --------------> B --------------------> C -------> Y
    //          ^ ----------------------------------------- D
    //                                2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{250.0, 0.0}});
    auto c = tg.addNd({{500.0, 0}});
    auto d = tg.addNd({{510.0, 5}});
    auto y = tg.addNd({{800.0, 1.0}});
    auto x = tg.addNd({{-300.0, 0.0}});

    auto xa = tg.addEdg(x, a, {{{-300.0, 0.0},  {0.0, 0.0}}});
    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {250.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{250.0, 0.0}, {500.0, 0}}});
    auto da = tg.addEdg(d, a, {{{510.0, 5.0}, {0.0, 0.0}}});
    auto cy = tg.addEdg(c, y, {{{500.0, 0.0}, {800.0, 1.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    xa->pl().addRoute(&l1, 0);
    xa->pl().addRoute(&l2, 0);

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    da->pl().addRoute(&l2, 0);
    cy->pl().addRoute(&l1, 0);
    cy->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;
    cfg.minSegLength = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();
    mc.removeNodeArtifacts();

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(tg.getNds()->size() == 2);

    a = x->getAdjList().front()->getOtherNd(x);
    d = y->getAdjList().front()->getOtherNd(y);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 2);
  }
}
