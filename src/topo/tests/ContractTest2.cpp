// Copyright 2016
// Author: Patrick Brosi


#include <cassert>
#include <string>

#include "shared/linegraph/LineGraph.h"
#include "topo/config/TopoConfig.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "topo/tests/ContractTest2.h"
#include "topo/tests/TopoTestUtil.h"
#include "util/geo/output/GeoGraphJsonOutput.h"

#define private public
#include "topo/mapconstructor/MapConstructor.h"

using util::approx;

// _____________________________________________________________________________
void ContractTest2::run() {
  // ___________________________________________________________________________
  {
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{0.0, -1.0}});
    auto c = tg.addNd({{2.0, 0.0}});
    auto d = tg.addNd({{-1.5, 0.0}});
    auto e = tg.addNd({{0.5, -0.5}});

    auto x = tg.addNd({{0, 200.0}});
    auto y = tg.addNd({{0, -200.0}});

    c->pl().addStop(shared::linegraph::Station("Thorndale", "Thorndale", *c->pl().getGeom()));

    auto ax = tg.addEdg(a, x, {{{0.0, 0.0}, {0.0, 200.0}}});
    auto yb = tg.addEdg(y, b, {{{0.0, -200.0}, {0.0, -1.0}}});

    auto ad = tg.addEdg(a, d, {{{0.0, 0.0}, {-1.5, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{0.0, -1.0}, {-1.5, 0.0}}});
    auto ce = tg.addEdg(c, e, {{{2.0, 0.0}, {0.5, -.5}}});
    auto dc = tg.addEdg(d, c, util::geo::PolyLine<double>({{-1.5, 0.0}, {0, 1}, {2.0, 0.0}}));
    auto ea = tg.addEdg(e, a, {{{0.5, -.5}, {0.0, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{0.5, -.5}, {0.0, -1.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "green");

    ax->pl().addLine(&l1, 0);
    ax->pl().addLine(&l2, 0);

    yb->pl().addLine(&l1, 0);
    yb->pl().addLine(&l2, 0);

    ad->pl().addLine(&l2, 0);
    bd->pl().addLine(&l1, 0);
    bd->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);
    dc->pl().addLine(&l1, 0);
    ea->pl().addLine(&l1, 0);
    ea->pl().addLine(&l2, 0);
    eb->pl().addLine(&l2, 0);


    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(a, d);
    mc.combineNodes(e, c);
    mc.combineNodes(d, b);
    mc.combineNodes(b, c);

  }
  // ___________________________________________________________________________
  {
    /*
     *
     *               e
     *              ^ ^
     *          ^  /   \
     *        1/  /     \ 1\
     *    ->     /  ->   \  v
     *    1     /    1    \    1->
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    e->pl().addStop(shared::linegraph::Station("1", "1", *e->pl().getGeom()));

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{25.0, 10.0}, {20.0, 0.0}}});
    auto ec = tg.addEdg(e, c, {{{25.0, 10.0}, {30.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "green");

    ab->pl().addLine(&l1, b);
    bc->pl().addLine(&l1, c);
    cd->pl().addLine(&l1, d);
    eb->pl().addLine(&l1, e);
    ec->pl().addLine(&l1, c);

    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l2, 0);
    eb->pl().addLine(&l2, 0);
    ec->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    // TODO: fix figure above and add tests
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, c
    //    1      2
    // a ---> b ---> c

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "green");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, c);
  }
}
