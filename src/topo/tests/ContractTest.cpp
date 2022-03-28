// Copyright 2016
// Author: Patrick Brosi

#include <cassert>
#include <string>

#include "shared/linegraph/LineGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/tests/ContractTest.h"
#include "topo/tests/TopoTestUtil.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/output/GeoGraphJsonOutput.h"

#define private public
#include "topo/mapconstructor/MapConstructor.h"

using util::approx;

// _____________________________________________________________________________
void ContractTest::run() {
  // ___________________________________________________________________________
  {
    /*
     *                f
     *                |
     *                |  2,3
     *                |
     *                v
     *                e
     *               ^ ^
     *              /   \
     *           2 /     \ 2,3->c
     *            /       \
     *     1,2   /    3->c \    1,2
     *  a -----> b ------> c ------> d
     *           \         /
     *            \1,3    / 1
     *             \     /
     *              \   /
     *               v v
     *                g
     *                |
     *                | 1,3
     *                v
     *                h
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto g = tg.addNd({{25.0, -10.0}});
    auto h = tg.addNd({{25.0, -20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    auto bg = tg.addEdg(b, g, {{{20.0, 0.0}, {25.0, -10.0}}});
    auto cg = tg.addEdg(c, g, {{{30.0, 0.0}, {25.0, -10.0}}});
    auto gh = tg.addEdg(g, h, {{{25.0, -10.0}, {25.0, -20.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);

    bc->pl().addLine(&l3, c);

    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);

    bg->pl().addLine(&l1, 0);
    bg->pl().addLine(&l3, 0);

    cg->pl().addLine(&l1, 0);

    gh->pl().addLine(&l1, 0);
    gh->pl().addLine(&l3, 0);

    ce->pl().addLine(&l2, 0);
    ce->pl().addLine(&l3, c);

    fe->pl().addLine(&l2, 0);
    fe->pl().addLine(&l3, 0);

    be->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, c);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    TEST(c->getAdjList().size(), ==, 4);

    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(f, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(g, h)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, g)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, g)->pl().hasLine(&l3));
    TEST(tg.getEdg(c, g)->pl().hasLine(&l1));

    TEST(tg.getEdg(e, c)->pl().hasLine(&l2));

    TEST(tg.getEdg(e, c)->pl().hasLine(&l3));
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2,3
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2,3
     *           /       \
     *    1,2   /    3    \    1,2
     * a -----> b ------> c ------> d
     *          \         /
     *           \1,3    / 1
     *            \     /
     *             \   /
     *              v v
     *               g
     *               |
     *               | 1,3
     *               v
     *               h
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto g = tg.addNd({{25.0, -10.0}});
    auto h = tg.addNd({{25.0, -20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    auto bg = tg.addEdg(b, g, {{{20.0, 0.0}, {25.0, -10.0}}});
    auto cg = tg.addEdg(c, g, {{{30.0, 0.0}, {25.0, -10.0}}});
    auto gh = tg.addEdg(g, h, {{{25.0, -10.0}, {25.0, -20.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);

    bc->pl().addLine(&l3, 0);

    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);

    bg->pl().addLine(&l1, 0);
    bg->pl().addLine(&l3, 0);

    cg->pl().addLine(&l1, 0);

    gh->pl().addLine(&l1, 0);
    gh->pl().addLine(&l3, 0);

    ce->pl().addLine(&l2, 0);
    ce->pl().addLine(&l3, 0);

    fe->pl().addLine(&l2, 0);
    fe->pl().addLine(&l3, 0);

    be->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l3, bc, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    TEST(c->getAdjList().size(), ==, 4);

    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(f, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(g, h)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, g)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, g)->pl().hasLine(&l3));
    TEST(tg.getEdg(c, g)->pl().hasLine(&l1));

    TEST(tg.getEdg(e, c)->pl().hasLine(&l2));

    TEST(tg.getEdg(e, c)->pl().hasLine(&l3));
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2,3
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2,3
     *           /       \
     *    1,2   /    3    \    1,2
     * a -----> b ------> c ------> d
     *          \         /
     *           \1,3    / 1
     *            \     /
     *             \   /
     *              v v
     *               g
     *               |
     *               | 1,3
     *               v
     *               h
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto g = tg.addNd({{25.0, -10.0}});
    auto h = tg.addNd({{25.0, -20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    auto bg = tg.addEdg(b, g, {{{20.0, 0.0}, {25.0, -10.0}}});
    auto cg = tg.addEdg(c, g, {{{30.0, 0.0}, {25.0, -10.0}}});
    auto gh = tg.addEdg(g, h, {{{25.0, -10.0}, {25.0, -20.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);

    bc->pl().addLine(&l3, 0);

    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);

    bg->pl().addLine(&l1, 0);
    bg->pl().addLine(&l3, 0);

    cg->pl().addLine(&l1, 0);

    gh->pl().addLine(&l1, 0);
    gh->pl().addLine(&l3, 0);

    ce->pl().addLine(&l2, 0);
    ce->pl().addLine(&l3, 0);

    fe->pl().addLine(&l2, 0);
    fe->pl().addLine(&l3, 0);

    be->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;
    //
    TEST(c->getAdjList().size(), ==, 4);

    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(f, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(g, h)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, g)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, g)->pl().hasLine(&l3));
    TEST(tg.getEdg(c, g)->pl().hasLine(&l1));

    TEST(tg.getEdg(e, c)->pl().hasLine(&l3));
    TEST(tg.getEdg(e, c)->pl().hasLine(&l2));
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2,3
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2,3
     *           /       \
     *    1,2   /    3    \    1,2
     * a -----> b ------> c ------> d
     *          \         /
     *           \1,3    / 1
     *            \     /
     *             \   /
     *              v v
     *               g
     *               |
     *               | 1,3
     *               v
     *               h
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto g = tg.addNd({{25.0, -10.0}});
    auto h = tg.addNd({{25.0, -20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    auto bg = tg.addEdg(b, g, {{{20.0, 0.0}, {25.0, -10.0}}});
    auto cg = tg.addEdg(c, g, {{{30.0, 0.0}, {25.0, -10.0}}});
    auto gh = tg.addEdg(g, h, {{{25.0, -10.0}, {25.0, -20.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);

    bc->pl().addLine(&l3, 0);

    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);

    bg->pl().addLine(&l1, 0);
    bg->pl().addLine(&l3, 0);

    cg->pl().addLine(&l1, 0);

    gh->pl().addLine(&l1, 0);
    gh->pl().addLine(&l3, 0);

    ce->pl().addLine(&l2, 0);
    ce->pl().addLine(&l3, 0);

    fe->pl().addLine(&l2, 0);
    fe->pl().addLine(&l3, 0);

    be->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, b);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;
    //
    TEST(b->getAdjList().size(), ==, 4);

    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(b, d)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(f, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(g, h)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(b, g)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(b, e)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(b, g)->pl().hasLine(&l3));
    TEST(tg.getEdg(b, g)->pl().hasLine(&l1));

    TEST(tg.getEdg(e, b)->pl().hasLine(&l3));
    TEST(tg.getEdg(e, b)->pl().hasLine(&l2));
  }

  // ___________________________________________________________________________
  {
    // node contraction of c, b
    //    1      1
    // a ---> b ---> c
    //        ^
    //        |  1
    //        d
    //           (d is a terminus for 1)

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});

    shared::linegraph::Line l1("1", "1", "red");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    // node contraction of c, b
    //    1      1
    // a ---> b ---> c
    //        ^
    //        |  1|
    //        d   v
    //           (d is a terminus for 1, but 1 doesn leave d)

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});

    shared::linegraph::Line l1("1", "1", "red");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, d);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    // node contraction of c, b
    //    1      1
    // a ---> b ---> c
    //        ^
    //        |  1|
    //        d   v
    //        |  (d is a terminus for 1, but 1 doesn leave d)
    //     2  |
    //        v
    //        e

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});
    auto e = tg.addNd({{50.0, -100.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});
    auto de = tg.addEdg(d, e, {{{50.0, -50.0}, {50.0, -100.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "green");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, d);
    de->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    // node contraction of c, b
    //    1      1
    // a ---> b ---> c
    //        ^
    //        |  1
    //        d
    //        |  (d is a terminus for 1)
    //     2  |
    //        v
    //        e

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});
    auto e = tg.addNd({{50.0, -100.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});
    auto de = tg.addEdg(d, e, {{{50.0, -50.0}, {50.0, -100.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "green");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, 0);
    de->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    // node contraction of c, b
    //    1      1
    // a ---> b ---> c
    //        ^
    //        |  1
    //        d
    //        |  (d is a terminus for 1 because of an exception)
    //     1  |
    //        v
    //        e

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});
    auto e = tg.addNd({{50.0, -100.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});
    auto de = tg.addEdg(d, e, {{{50.0, -50.0}, {50.0, -100.0}}});

    shared::linegraph::Line l1("1", "1", "red");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, 0);
    de->pl().addLine(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    d->pl().addConnExc(&l1, de, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    /*
     * node contraction of c, b
     *    1      1
     * a ---> b ---> c
     *        ^    /
     *        | 1 / 1
     *        d -/
     *               (d is a terminus for 1 because of an exception)
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, 0);
    dc->pl().addLine(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    d->pl().addConnExc(&l1, dc, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    /*
     * node contraction of c, b
     *    1      1
     * a ---> b ---> c
     *        ^    /
     *        | 1 / 1
     *        d -/
     *               (d is a terminus for 1 because of an exception)
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, 0);
    dc->pl().addLine(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    d->pl().addConnExc(&l1, dc, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(c, d);

    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(d, b)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    /*
     * node contraction of c, b
     *    1      1
     * a ---> b ---> c
     *        ^    /
     *        | 1 / 2
     *        d -/
     *
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "green");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, 0);
    dc->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(c, d);

    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(d, b)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    /*
     * node contraction of c, b
     *    1      1
     * a ---> b ---> c
     *        ^    /
     *        | 1 / 2
     *        d -/
     *
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "green");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l1, 0);
    bd->pl().addLine(&l1, 0);
    dc->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);
    b->pl().addConnExc(&l1, ab, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(c, d);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(d, b)->pl().getLines().size(), ==, 1);
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);
    fe->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    TEST(tg.getNds().size(), ==, 5);
    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, e)->pl().getLines().begin()->line, ==, &l2);
    TEST(tg.getEdg(c, e)->pl().getLines().begin()->direction, ==, 0);
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);
    fe->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, b);

    TEST(tg.getNds().size(), ==, 5);
    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(b, d)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(b, e)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(b, e)->pl().getLines().begin()->line, ==, &l2);
    TEST(tg.getEdg(b, e)->pl().getLines().begin()->direction, ==, 0);
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);
    fe->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 5);
    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(b, e)->pl().getLines().size(), ==, 2);
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);
    fe->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 5);
    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(b, c)->pl().getLines().size(), ==, 2);
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);
    fe->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, b);

    TEST(tg.getNds().size(), ==, 5);
    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(f, b)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(f, b)->pl().getLines().begin()->line, ==, &l2);
    TEST(tg.getEdg(f, b)->pl().getLines().begin()->direction, ==, 0);

    TEST(tg.getEdg(b, c)->pl().getLines().size(), ==, 2);
  }

  // ___________________________________________________________________________
  {
    /*
     *               f
     *               |
     *               |  2
     *               |
     *               v
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});
    auto f = tg.addNd({{25.0, 20.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});
    auto fe = tg.addEdg(f, e, {{{25.0, 20.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);
    fe->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    TEST(tg.getNds().size(), ==, 5);
    TEST(tg.getEdg(a, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(f, e)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(f, e)->pl().getLines().begin()->line, ==, &l2);
    TEST(tg.getEdg(f, e)->pl().getLines().begin()->direction, ==, 0);

    TEST(tg.getEdg(e, c)->pl().getLines().size(), ==, 2);
  }

  // ==========================================================================

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));

    TEST(d->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(d->getAdjList().front()->pl().lineOcc(&l2).direction, ==,
         d->getAdjList().front()->getOtherNd(d));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    TEST(tg.getNds().size(), ==, 4);
    auto test = tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d));
    TEST(test->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(test->pl().hasLine(&l1));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2->   \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, c);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    TEST(tg.getNds().size(), ==, 4);

    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l3));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .lineOcc(&l2)
             .direction,
         ==, d->getAdjList().front()->getOtherNd(d));
    TEST(tg.getEdg(d, d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .lineOcc(&l2)
             .direction,
         ==, d->getAdjList().front()->getOtherNd(d));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \   ^
     *          3 /     \ 2|
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, e);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l3));
  }

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
  }

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getLines().size() == 0)
          e = edg->getOtherNd(c);
        else
          a = edg->getOtherNd(c);
      }
    }

    TEST(tg.getNds().size(), ==, 4);

    TEST(tg.getEdg(a, c));
    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, d));
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, e));
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 1);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getLines().size() == 0)
          e = edg->getOtherNd(c);
        else
          a = edg->getOtherNd(c);
      }
    }

    TEST(tg.getNds().size(), ==, 4);

    TEST(tg.getEdg(a, c));
    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, d));
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, e));
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(c, e)->pl().hasLine(&l2));
  }

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2->   \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, c);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getLines().size() == 0)
          e = edg->getOtherNd(c);
        else
          a = edg->getOtherNd(c);
      }
    }

    TEST(tg.getNds().size(), ==, 4);
    // all 3 lines in the triangle are lost
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getLines().size() == 0)
          e = edg->getOtherNd(c);
        else
          a = edg->getOtherNd(c);
      }
    }

    TEST(tg.getNds().size(), ==, 4);
    // all 3 lines in the triangle are lost
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getLines().size() == 0)
          e = edg->getOtherNd(c);
        else
          a = edg->getOtherNd(c);
      }
    }

    // bc will now be contracted and be and ce will be folded into one edge
    // with no lines on it, because 2 and 3 are lost on it
    // neither 1 nor 2 can continue from ab to cd because 1 cannot continue
    // through the triangle, and 2 is restricted by the exception from cd to bc

    TEST(tg.getNds().size(), ==, 4);

    TEST(tg.getEdg(a, c));
    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, d));
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, e));
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // bc will now be contracted and be and ce will be folded into one edge
    // with no lines on it, because 2 and 1 are lost on it

    TEST(c->getAdjList().size(), ==, 3);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getLines().size() == 0)
          e = edg->getOtherNd(c);
        else
          a = edg->getOtherNd(c);
      }
    }

    TEST(tg.getNds().size(), ==, 4);

    TEST(tg.getEdg(a, c));
    TEST(tg.getEdg(a, c)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, d));
    TEST(tg.getEdg(c, d)->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(c, e));
    TEST(tg.getEdg(c, e)->pl().getLines().size(), ==, 2);
  }

  // ===========================================================================

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    // line 1 on bc is deleted because it cannot leave bc
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
    TEST(d->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(d->getAdjList().front()->pl().lineOcc(&l2).direction, ==,
         d->getAdjList().front()->getOtherNd(d));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    // line 2 is deleted because it is lost
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2->   \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, c);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);

    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l3));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    // 1 on ec was deleted, 2 on be was deleted
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
  }

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    // line 1 on bc is deleted because it cannot leave bc
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));

    TEST(d->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(d->getAdjList().front()->pl().lineOcc(&l2).direction, ==,
         d->getAdjList().front()->getOtherNd(d));
  }

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    // line 2 is deleted because it is lost
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2->   \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, c);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);

    // line 3 and line 2-> on bc are dead ends (because ec is contracted), line
    // 2 on ec is also a dead end
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
  }

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    // all 3 lines in the triangle are lost
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);

    // 3 on be and 2 on bc are lost, 2 on ce is contracted but continues outside
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l3));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
  }

  // ===========================================================================

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));

    TEST(d->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(d->getAdjList().front()->pl().lineOcc(&l2).direction, ==,
         d->getAdjList().front()->getOtherNd(d));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    // line 2 is deleted because it is lost
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2->   \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, c);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);

    // line 3 and line 2-> on bc are dead ends (because ec is contracted), line
    // 2 on ec is also a dead end
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    // all 3 lines in the triangle are lost
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);

    // 3 on be and 2 on bc are lost, 2 on ce is contracted but continues outside
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l3));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    // 1 on ec was deleted, 2 on be was deleted
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
  }

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    // line 2 is deleted because it is lost
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2->   \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, c);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);

    // line 3 and line 2-> on bc are dead ends (because ec is contracted), line
    // 2 on ec is also a dead end
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    // all 3 lines in the triangle are lost
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);

    // 3 on be and 2 on bc are lost, 2 on ce is contracted but continues outside
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l3));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
  }

  // ===========================================================================

  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->direction,
         ==, 0);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 2
     *           /       \
     *    1,2   /   1     \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l1, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    TEST(tg.getNds().size(), ==, 4);

    // 2 on be and 2 on ce are dead ends
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->direction,
         ==, 0);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2->   \    1,<-2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, c);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, c);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->direction,
         ==, 0);

    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(d->getAdjList().front()->pl().lineOcc(&l2).direction, ==,
         d->getAdjList().front()->getOtherNd(d));
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          3 /     \ 2
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l3, 0);
    ce->pl().addLine(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->direction,
         ==, 0);

    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l2));
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .hasLine(&l1));

    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .lineOcc(&l2)
             .direction,
         ==, 0);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .lineOcc(&l1)
             .direction,
         ==, 0);
  }
  // ___________________________________________________________________________
  {
    /*
     *               e
     *              ^ ^
     *             /   \
     *          2 /     \ 1
     *           /       \
     *    1,2   /   2     \    1,2
     * a -----> b ------> c ------> d
     */

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{10.0, 0.0}});
    auto b = tg.addNd({{20.0, 0.0}});
    auto c = tg.addNd({{30.0, 0.0}});
    auto d = tg.addNd({{40.0, 0.0}});
    auto e = tg.addNd({{25.0, 10.0}});

    auto ab = tg.addEdg(a, b, {{{10.0, 0.0}, {20.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{20.0, 0.0}, {30.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{30.0, 0.0}, {40.0, 0.0}}});
    auto be = tg.addEdg(b, e, {{{20.0, 0.0}, {25.0, 10.0}}});
    auto ce = tg.addEdg(c, e, {{{30.0, 0.0}, {25.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);
    be->pl().addLine(&l2, 0);
    ce->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    TEST(tg.getNds().size(), ==, 4);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .size(),
         ==, 1);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->line,
         ==, &l2);
    TEST(tg.getEdg(a->getAdjList().front()->getOtherNd(a),
                   d->getAdjList().front()->getOtherNd(d))
             ->pl()
             .getLines()
             .begin()
             ->direction,
         ==, 0);
  }

  // ___________________________________________________________________________
  {
    //     1        2          1
    // a -----> b ------> c ------> d
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{100.0, 0.0}});
    auto b = tg.addNd({{200.0, 0.0}});
    auto c = tg.addNd({{300.0, 0.0}});
    auto d = tg.addNd({{400.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{100.0, 0.0}, {200.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{200.0, 0.0}, {300.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{300.0, 0.0}, {400.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    TEST(tg.getNds().size(), ==, 3);
    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 1);

    TEST(a->getAdjList().front()->pl().getLines().begin()->line, ==, &l1);
    TEST(d->getAdjList().front()->pl().getLines().begin()->line, ==, &l1);

    TEST(a->getAdjList().front()->pl().getLines().begin()->direction, ==, 0);
    TEST(d->getAdjList().front()->pl().getLines().begin()->direction, ==, 0);
  }

  // ___________________________________________________________________________
  {
    //    1         2          1
    // a -----> b ------> c ------> d
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{100.0, 0.0}});
    auto b = tg.addNd({{200.0, 0.0}});
    auto c = tg.addNd({{300.0, 0.0}});
    auto d = tg.addNd({{400.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{100.0, 0.0}, {200.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{200.0, 0.0}, {300.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{300.0, 0.0}, {400.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    bc->pl().addLine(&l2, 0);
    cd->pl().addLine(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, b);

    TEST(tg.getNds().size(), ==, 3);
    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 1);

    TEST(a->getAdjList().front()->pl().getLines().begin()->line, ==, &l1);
    TEST(d->getAdjList().front()->pl().getLines().begin()->line, ==, &l1);

    TEST(a->getAdjList().front()->pl().getLines().begin()->direction, ==, 0);
    TEST(d->getAdjList().front()->pl().getLines().begin()->direction, ==, 0);
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c ---> a ---> b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ca->pl().addLine(&l1, 0);
    ca->pl().addLine(&l2, 0);
    ca->pl().addLine(&l3, 0);

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    ab->pl().addLine(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ca, ab));

    mc.combineEdges(ca, ab, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c <--- a ---> b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});
    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    ab->pl().addLine(&l3, 0);

    ac->pl().addLine(&l1, 0);
    ac->pl().addLine(&l2, 0);
    ac->pl().addLine(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ab, ac));

    mc.combineEdges(ab, ac, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c <--- a <--- b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ba->pl().addLine(&l1, 0);
    ba->pl().addLine(&l2, 0);
    ba->pl().addLine(&l3, 0);

    ac->pl().addLine(&l1, 0);
    ac->pl().addLine(&l2, 0);
    ac->pl().addLine(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ac, ba));

    mc.combineEdges(ac, ba, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c ---> a <--- b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ba->pl().addLine(&l1, 0);
    ba->pl().addLine(&l2, 0);
    ba->pl().addLine(&l3, 0);

    ca->pl().addLine(&l1, 0);
    ca->pl().addLine(&l2, 0);
    ca->pl().addLine(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ca, ba));

    mc.combineEdges(ca, ba, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // with line directions and without connection exclusions
    // c ---> a ---> b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ca->pl().addLine(&l1, 0);
    ca->pl().addLine(&l2, a);
    ca->pl().addLine(&l3, c);

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, b);
    ab->pl().addLine(&l3, a);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ca, ab));

    mc.combineEdges(ca, ab, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, b);
      if (r.line == &l3) TEST(r.direction, ==, c);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c <--- a ---> b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});
    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, b);
    ab->pl().addLine(&l3, a);

    ac->pl().addLine(&l1, 0);
    ac->pl().addLine(&l2, a);
    ac->pl().addLine(&l3, c);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ac, ab));

    mc.combineEdges(ab, ac, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, b);
      if (r.line == &l3) TEST(r.direction, ==, c);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c <--- a <--- b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ba->pl().addLine(&l1, 0);
    ba->pl().addLine(&l2, b);
    ba->pl().addLine(&l3, a);

    ac->pl().addLine(&l1, 0);
    ac->pl().addLine(&l2, a);
    ac->pl().addLine(&l3, c);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ba, ac));

    mc.combineEdges(ac, ba, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, b);
      if (r.line == &l3) TEST(r.direction, ==, c);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c ---> a <--- b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ba->pl().addLine(&l1, 0);
    ba->pl().addLine(&l2, b);
    ba->pl().addLine(&l3, a);

    ca->pl().addLine(&l1, 0);
    ca->pl().addLine(&l2, a);
    ca->pl().addLine(&l3, c);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ba, ca));

    mc.combineEdges(ca, ba, a);

    TEST(tg.getNds().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(b->getAdjList().front(), ==, c->getAdjList().front());
    TEST(b->getAdjList().front()->pl().getLines().size(), ==, 3);

    // check geometry orientation
    TEST(b->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(b->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, b);
      if (r.line == &l3) TEST(r.direction, ==, c);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c ---> a ---> b
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    TEST(tg.getNds().size(), ==, 3);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 1);
    TEST(c->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ca->pl().addLine(&l1, 0);
    ca->pl().addLine(&l2, 0);
    ca->pl().addLine(&l3, 0);

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    ab->pl().addLine(&l3, 0);

    a->pl().addConnExc(&l1, ca, ab);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(!mc.lineEq(ca, ab));
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c ---> a ---> b ---> d
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    TEST(tg.getNds().size(), ==, 4);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ca->pl().addLine(&l1, 0);
    ca->pl().addLine(&l2, 0);
    ca->pl().addLine(&l3, 0);

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    ab->pl().addLine(&l3, 0);

    bd->pl().addLine(&l3, 0);
    bd->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l3, ab, bd);
    b->pl().addConnExc(&l2, bd, ab);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ca, ab));

    mc.combineEdges(ca, ab, a);

    TEST(tg.getNds().size(), ==, 3);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 3);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);

    // check geometry orientation
    TEST(c->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(c->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c <--- a ---> b ---> d
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    TEST(tg.getNds().size(), ==, 4);

    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ac->pl().addLine(&l1, 0);
    ac->pl().addLine(&l2, 0);
    ac->pl().addLine(&l3, 0);

    ab->pl().addLine(&l1, 0);
    ab->pl().addLine(&l2, 0);
    ab->pl().addLine(&l3, 0);

    bd->pl().addLine(&l3, 0);
    bd->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l3, ab, bd);
    b->pl().addConnExc(&l2, bd, ab);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ac, ab));

    mc.combineEdges(ac, ab, a);

    TEST(tg.getNds().size(), ==, 3);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 3);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);

    // check geometry orientation
    TEST(c->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(c->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c ---> a <--- b ---> d
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    TEST(tg.getNds().size(), ==, 4);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ca->pl().addLine(&l1, 0);
    ca->pl().addLine(&l2, 0);
    ca->pl().addLine(&l3, 0);

    ba->pl().addLine(&l1, 0);
    ba->pl().addLine(&l2, 0);
    ba->pl().addLine(&l3, 0);

    bd->pl().addLine(&l3, 0);
    bd->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l3, ba, bd);
    b->pl().addConnExc(&l2, bd, ba);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ba, ba));

    mc.combineEdges(ca, ba, a);

    TEST(tg.getNds().size(), ==, 3);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 3);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);

    // check geometry orientation
    TEST(c->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(c->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c <--- a <--- b ---> d
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    TEST(tg.getNds().size(), ==, 4);

    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});
    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    TEST(a->getAdjList().size(), ==, 2);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ac->pl().addLine(&l1, 0);
    ac->pl().addLine(&l2, 0);
    ac->pl().addLine(&l3, 0);

    ba->pl().addLine(&l1, 0);
    ba->pl().addLine(&l2, 0);
    ba->pl().addLine(&l3, 0);

    bd->pl().addLine(&l3, 0);
    bd->pl().addLine(&l2, 0);

    b->pl().addConnExc(&l3, ba, bd);
    b->pl().addConnExc(&l2, bd, ba);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    TEST(mc.lineEq(ba, ba));

    mc.combineEdges(ac, ba, a);

    TEST(tg.getNds().size(), ==, 3);
    TEST(b->getAdjList().size(), ==, 2);
    TEST(c->getAdjList().size(), ==, 1);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 3);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 2);

    // check geometry orientation
    TEST(c->getAdjList().front()->pl().getGeom()->front().getX(), ==,
         approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    TEST(c->getAdjList().front()->pl().getGeom()->back().getX(), ==,
         approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, d
    //        c
    //     2,3^
    //   1    | 1,2,3  1,2,3
    // a ---> b ---> d ---> e

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});

    TEST(tg.getNds().size(), ==, 5);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});

    TEST(a->getAdjList().size(), ==, 1);
    TEST(b->getAdjList().size(), ==, 3);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 2);
    TEST(e->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);

    bc->pl().addLine(&l2, 0);
    bc->pl().addLine(&l3, 0);

    bd->pl().addLine(&l1, 0);
    bd->pl().addLine(&l2, 0);
    bd->pl().addLine(&l3, 0);

    de->pl().addLine(&l1, 0);
    de->pl().addLine(&l2, 0);
    de->pl().addLine(&l3, 0);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    // should now be
    //        c
    //     2,3^
    //    1   |  1,2,3
    // a ---> d ---> e

    TEST(tg.getNds().size(), ==, 4);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 3);
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, d
    //        c
    //     2,3^
    //   1    | 1,2,3  1,2,3
    // a ---> b ---> d ---> e

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});

    TEST(tg.getNds().size(), ==, 5);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});

    TEST(a->getAdjList().size(), ==, 1);
    TEST(b->getAdjList().size(), ==, 3);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 2);
    TEST(e->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, b);

    bc->pl().addLine(&l2, c);
    bc->pl().addLine(&l3, c);

    bd->pl().addLine(&l1, d);
    bd->pl().addLine(&l2, b);
    bd->pl().addLine(&l3, b);

    de->pl().addLine(&l1, e);
    de->pl().addLine(&l2, d);
    de->pl().addLine(&l3, d);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    // should now be
    //        c
    //     2,3^
    //    1   |  1,2,3
    // a ---> d ---> e

    TEST(tg.getNds().size(), ==, 4);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 3);

    for (auto r : e->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, e);
      if (r.line == &l2) TEST(r.direction, ==, d);
      if (r.line == &l3) TEST(r.direction, ==, d);
    }

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      if (r.line == &l2) TEST(r.direction, ==, c);
      if (r.line == &l3) TEST(r.direction, ==, c);
    }

    for (auto r : a->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, d);
    }
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, d
    //        c
    //     2,3^
    //   1    | 1,2,3  1,2,3
    // a ---> b ---> d ---> e

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});

    TEST(tg.getNds().size(), ==, 5);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});

    TEST(a->getAdjList().size(), ==, 1);
    TEST(b->getAdjList().size(), ==, 3);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 2);
    TEST(e->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);

    bc->pl().addLine(&l2, 0);
    bc->pl().addLine(&l3, 0);

    bd->pl().addLine(&l1, 0);
    bd->pl().addLine(&l2, 0);
    bd->pl().addLine(&l3, 0);

    de->pl().addLine(&l1, 0);
    de->pl().addLine(&l2, 0);
    de->pl().addLine(&l3, 0);

    // no connection of line 2 from cb -> bd
    b->pl().addConnExc(&l2, bc, bd);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    // should now be
    //        c
    //     2,3^
    //    1   |  1,2,3
    // a ---> d ---> e

    TEST(tg.getNds().size(), ==, 4);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 3);

    for (auto r : e->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, 0);
      if (r.line == &l3) TEST(r.direction, ==, 0);
    }

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      if (r.line == &l2) TEST(r.direction, ==, 0);
      if (r.line == &l3) TEST(r.direction, ==, 0);
    }

    for (auto r : a->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, d
    //        c
    //     2,3^
    //   1    | 1,2,3  1,3
    // a ---> b ---> d ---> e
    //               |  2
    //               -----> f

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});
    auto f = tg.addNd({{150.0, -50.0}});

    TEST(tg.getNds().size(), ==, 6);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});
    auto df = tg.addEdg(d, f, {{{100.0, 0.0}, {150.0, -50.0}}});

    TEST(a->getAdjList().size(), ==, 1);
    TEST(b->getAdjList().size(), ==, 3);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 3);
    TEST(e->getAdjList().size(), ==, 1);
    TEST(f->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);

    bc->pl().addLine(&l2, 0);
    bc->pl().addLine(&l3, 0);

    bd->pl().addLine(&l1, 0);
    bd->pl().addLine(&l2, 0);
    bd->pl().addLine(&l3, 0);

    de->pl().addLine(&l1, 0);
    de->pl().addLine(&l3, 0);

    df->pl().addLine(&l2, 0);

    // no connection of line 3 from cb -> bd
    b->pl().addConnExc(&l3, bc, bd);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    // should now be
    //        c
    //     2,3^
    //   1    | 1,3
    // a ---> d ---> e
    //        |  2
    //        -----> f

    TEST(tg.getNds().size(), ==, 5);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(f->getAdjList().front()->pl().getLines().size(), ==, 1);

    for (auto r : e->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l3) TEST(r.direction, ==, 0);
    }

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      if (r.line == &l2) TEST(r.direction, ==, 0);
      if (r.line == &l3) TEST(r.direction, ==, 0);
    }

    for (auto r : a->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
    }

    for (auto r : f->getAdjList().front()->pl().getLines()) {
      if (r.line == &l2) TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    // node contraction of d, b
    //        c
    //     2,3^
    //   1    | 1,2,3  1,3
    // a ---> b ---> d ---> e
    //               |  2
    //               -----> f

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});
    auto f = tg.addNd({{150.0, -50.0}});

    TEST(tg.getNds().size(), ==, 6);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});
    auto df = tg.addEdg(d, f, {{{100.0, 0.0}, {150.0, -50.0}}});

    TEST(a->getAdjList().size(), ==, 1);
    TEST(b->getAdjList().size(), ==, 3);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(d->getAdjList().size(), ==, 3);
    TEST(e->getAdjList().size(), ==, 1);
    TEST(f->getAdjList().size(), ==, 1);

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");

    ab->pl().addLine(&l1, 0);

    bc->pl().addLine(&l2, 0);
    bc->pl().addLine(&l3, 0);

    bd->pl().addLine(&l1, 0);
    bd->pl().addLine(&l2, 0);
    bd->pl().addLine(&l3, 0);

    de->pl().addLine(&l1, 0);
    de->pl().addLine(&l3, 0);

    df->pl().addLine(&l2, 0);

    // no connection of line 3 from cb -> bd
    b->pl().addConnExc(&l3, bc, bd);

    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(d, b);

    // should now be
    //        c
    //     2,3^
    //   1    | 1,3
    // a ---> b ---> e
    //        |  2
    //        -----> f

    TEST(tg.getNds().size(), ==, 5);
    TEST(c->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().size(), ==, 1);
    TEST(a->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);
    TEST(f->getAdjList().front()->pl().getLines().size(), ==, 1);

    for (auto r : e->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l3) TEST(r.direction, ==, 0);
    }

    for (auto r : c->getAdjList().front()->pl().getLines()) {
      if (r.line == &l2) TEST(r.direction, ==, 0);
      if (r.line == &l3) TEST(r.direction, ==, 0);
    }

    for (auto r : a->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
    }

    for (auto r : f->getAdjList().front()->pl().getLines()) {
      if (r.line == &l2) TEST(r.direction, ==, 0);
    }

    //
  }

  // ___________________________________________________________________________
  {
    //         1
    //     a -----> b
    //  1  |        ^
    //     c -------|
    //          2

    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {100.0, 10.0}}});
    auto ac = tg.addEdg(a, c, {{{0.0, 10.0}, {0.0, 5.0}}});
    auto cb = tg.addEdg(c, b, {{{0.0, 5.0}, {100.0, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    ac->pl().addLine(&l1, 0);
    cb->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.removeEdgeArtifacts();

    // this should've done nothing - it may look like we can contract a and c
    // and merge edges ab and cb into one, but in general, we cannot be sure
    // that this is even possible - their geometries may differ vastly. This
    // should be handled by the edge merging mechanism, and there should be
    // checks which prevent such a contraction
  }
}
