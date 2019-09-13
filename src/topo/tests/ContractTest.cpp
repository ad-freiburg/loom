// Copyright 2016
// Author: Patrick Brosi


#include <cassert>
#include <string>

#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "topo/tests/ContractTest.h"
#include "topo/tests/TopoTestUtil.h"
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);

    bc->pl().addRoute(&l3, c);

    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);

    bg->pl().addRoute(&l1, 0);
    bg->pl().addRoute(&l3, 0);

    cg->pl().addRoute(&l1, 0);

    gh->pl().addRoute(&l1, 0);
    gh->pl().addRoute(&l3, 0);

    ce->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l3, c);

    fe->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l3, 0);

    be->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, c);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(c->getAdjList().size() == 4);

    assert(tg.getEdg(a,c)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,d)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(f,e)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(g,h)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,g)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,e)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c,g)->pl().hasRoute(&l3));
    assert(tg.getEdg(c,g)->pl().hasRoute(&l1));

    assert(tg.getEdg(e,c)->pl().hasRoute(&l2));

    assert(tg.getEdg(e,c)->pl().hasRoute(&l3));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);

    bc->pl().addRoute(&l3, 0);

    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);

    bg->pl().addRoute(&l1, 0);
    bg->pl().addRoute(&l3, 0);

    cg->pl().addRoute(&l1, 0);

    gh->pl().addRoute(&l1, 0);
    gh->pl().addRoute(&l3, 0);

    ce->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l3, 0);

    fe->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l3, 0);

    be->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l3, bc, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(c->getAdjList().size() == 4);

    assert(tg.getEdg(a,c)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,d)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(f,e)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(g,h)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,g)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,e)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c,g)->pl().hasRoute(&l3));
    assert(tg.getEdg(c,g)->pl().hasRoute(&l1));

    assert(tg.getEdg(e,c)->pl().hasRoute(&l2));

    assert(tg.getEdg(e,c)->pl().hasRoute(&l3));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);

    bc->pl().addRoute(&l3, 0);

    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);

    bg->pl().addRoute(&l1, 0);
    bg->pl().addRoute(&l3, 0);

    cg->pl().addRoute(&l1, 0);

    gh->pl().addRoute(&l1, 0);
    gh->pl().addRoute(&l3, 0);

    ce->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l3, 0);

    fe->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l3, 0);

    be->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;
    //
    assert(c->getAdjList().size() == 4);

    assert(tg.getEdg(a,c)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,d)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(f,e)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(g,h)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,g)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c,e)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c,g)->pl().hasRoute(&l3));
    assert(tg.getEdg(c,g)->pl().hasRoute(&l1));

    assert(tg.getEdg(e,c)->pl().hasRoute(&l3));
    assert(tg.getEdg(e,c)->pl().hasRoute(&l2));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);

    bc->pl().addRoute(&l3, 0);

    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);

    bg->pl().addRoute(&l1, 0);
    bg->pl().addRoute(&l3, 0);

    cg->pl().addRoute(&l1, 0);

    gh->pl().addRoute(&l1, 0);
    gh->pl().addRoute(&l3, 0);

    ce->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l3, 0);

    fe->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l3, 0);

    be->pl().addRoute(&l2, 0);


    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, b);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;
    //
    assert(b->getAdjList().size() == 4);

    assert(tg.getEdg(a,b)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(b,d)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(f,e)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(g,h)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(b,g)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(b,e)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(b,g)->pl().hasRoute(&l3));
    assert(tg.getEdg(b,g)->pl().hasRoute(&l1));

    assert(tg.getEdg(e,b)->pl().hasRoute(&l3));
    assert(tg.getEdg(e,b)->pl().hasRoute(&l2));
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});

    transitmapper::graph::Route l1("1", "1", "red");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    assert(tg.getEdg(a, d)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});

    transitmapper::graph::Route l1("1", "1", "red");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, d);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    assert(tg.getEdg(a, d)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});
    auto e = tg.addNd({{50.0, -100.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});
    auto de = tg.addEdg(d, e, {{{50.0, -50.0}, {50.0, -100.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, d);
    de->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    assert(tg.getEdg(a, d)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});
    auto e = tg.addNd({{50.0, -100.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});
    auto de = tg.addEdg(d, e, {{{50.0, -50.0}, {50.0, -100.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    assert(tg.getEdg(a, d)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{100.0, 0.0}});
    auto d = tg.addNd({{50.0, -50.0}});
    auto e = tg.addNd({{50.0, -100.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {50.0, -50.0}}});
    auto de = tg.addEdg(d, e, {{{50.0, -50.0}, {50.0, -100.0}}});

    transitmapper::graph::Route l1("1", "1", "red");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    d->pl().addConnExc(&l1, de, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    assert(tg.getEdg(a, d)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, 0);
    dc->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    d->pl().addConnExc(&l1, dc, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    assert(tg.getEdg(a, d)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, 0);
    dc->pl().addRoute(&l1, 0);

    b->pl().addConnExc(&l1, ab, bc);
    d->pl().addConnExc(&l1, dc, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(c, d);

    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(d, b)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, 0);
    dc->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(c, d);


    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(d, b)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{5.0, 0.0}});
    auto c = tg.addNd({{10.0, 0.0}});
    auto d = tg.addNd({{5.0, -5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {5.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{5.0, 0.0}, {10.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{5.0, 0.0}, {5.0, -5.0}}});
    auto dc = tg.addEdg(d, c, {{{5.0, -5.0}, {10.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l1, 0);
    dc->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);
    b->pl().addConnExc(&l1, ab, bd);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(c, d);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(d, b)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    assert(tg.getNds()->size() == 5);
    assert(tg.getEdg(a, c)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(c, e)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, e)->pl().getRoutes().begin()->route == &l2);
    assert(tg.getEdg(c, e)->pl().getRoutes().begin()->direction == 0);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, b);

    assert(tg.getNds()->size() == 5);
    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(b, d)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(b, e)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(b, e)->pl().getRoutes().begin()->route == &l2);
    assert(tg.getEdg(b, e)->pl().getRoutes().begin()->direction == 0);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 5);
    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(b, e)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    assert(tg.getNds()->size() == 5);
    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(b, c)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, b);

    assert(tg.getNds()->size() == 5);
    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(f, b)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(f, b)->pl().getRoutes().begin()->route == &l2);
    assert(tg.getEdg(f, b)->pl().getRoutes().begin()->direction == 0);

    assert(tg.getEdg(b, c)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);
    fe->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, ce, be);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    assert(tg.getNds()->size() == 5);
    assert(tg.getEdg(a, e)->pl().getRoutes().size() == 2);
    assert(tg.getEdg(f, e)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(f, e)->pl().getRoutes().begin()->route == &l2);
    assert(tg.getEdg(f, e)->pl().getRoutes().begin()->direction == 0);

    assert(tg.getEdg(e, c)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;


    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));

    assert(d->getAdjList().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(d->getAdjList().front()->pl().getRouteOcc(&l2).direction == d->getAdjList().front()->getOtherNd(d));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l1);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, c);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(tg.getNds()->size() == 4);

    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l3));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRouteOcc(&l2).direction == d->getAdjList().front()->getOtherNd(d));
    assert(tg.getEdg(d, d->getAdjList().front()->getOtherNd(d))->pl().getRouteOcc(&l2).direction == d->getAdjList().front()->getOtherNd(d));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, e);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);


    assert(tg.getNds()->size() == 4);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l3));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;


    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);


    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getRoutes().size() == 0) e = edg->getOtherNd(c);
        else a = edg->getOtherNd(c);
      }
    }

    assert(tg.getNds()->size() == 4);

    assert(tg.getEdg(a, c));
    assert(tg.getEdg(a, c)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, d));
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, e));
    assert(tg.getEdg(c, e)->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);


    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getRoutes().size() == 0) e = edg->getOtherNd(c);
        else a = edg->getOtherNd(c);
      }
    }

    assert(tg.getNds()->size() == 4);

    assert(tg.getEdg(a, c));
    assert(tg.getEdg(a, c)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, d));
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, e));
    assert(tg.getEdg(c, e)->pl().getRoutes().size() == 1);
    assert(tg.getEdg(c, e)->pl().hasRoute(&l2));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, c);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getRoutes().size() == 0) e = edg->getOtherNd(c);
        else a = edg->getOtherNd(c);
      }
    }

    assert(tg.getNds()->size() == 4);
    // all 3 lines in the triangle are lost
    assert(tg.getEdg(c, e)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getRoutes().size() == 0) e = edg->getOtherNd(c);
        else a = edg->getOtherNd(c);
      }
    }

    assert(tg.getNds()->size() == 4);
    // all 3 lines in the triangle are lost
    assert(tg.getEdg(c, e)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);


    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getRoutes().size() == 0) e = edg->getOtherNd(c);
        else a = edg->getOtherNd(c);
      }
    }

    // bc will now be contracted and be and ce will be folded into one edge
    // with no lines on it, because 2 and 3 are lost on it
    // neither 1 nor 2 can continue from ab to cd because 1 cannot continue
    // through the triangle, and 2 is restricted by the exception from cd to bc

    assert(tg.getNds()->size() == 4);

    assert(tg.getEdg(a, c));
    assert(tg.getEdg(a, c)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, d));
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, e));
    assert(tg.getEdg(c, e)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    // bc will now be contracted and be and ce will be folded into one edge
    // with no lines on it, because 2 and 1 are lost on it


    assert(c->getAdjList().size() == 3);

    // c will be the kept node
    for (auto edg : c->getAdjList()) {
      if (edg->getOtherNd(c) != d) {
        if (edg->pl().getRoutes().size() == 0) e = edg->getOtherNd(c);
        else a = edg->getOtherNd(c);
      }
    }

    assert(tg.getNds()->size() == 4);

    assert(tg.getEdg(a, c));
    assert(tg.getEdg(a, c)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, d));
    assert(tg.getEdg(c, d)->pl().getRoutes().size() == 2);

    assert(tg.getEdg(c, e));
    assert(tg.getEdg(c, e)->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;


    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    // line 1 on bc is deleted because it cannot leave bc
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
    assert(d->getAdjList().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(d->getAdjList().front()->pl().getRouteOcc(&l2).direction == d->getAdjList().front()->getOtherNd(d));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    // line 2 is deleted because it is lost
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, c);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);


    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);


    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);

    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l3));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    // 1 on ec was deleted, 2 on be was deleted
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    // line 1 on bc is deleted because it cannot leave bc
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));

    assert(d->getAdjList().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(d->getAdjList().front()->pl().getRouteOcc(&l2).direction == d->getAdjList().front()->getOtherNd(d));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    // line 2 is deleted because it is lost
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, c);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);

    // line 3 and line 2-> on bc are dead ends (because ec is contracted), line 2 on ec is also a dead end
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
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
    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    // all 3 lines in the triangle are lost
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);

    // 3 on be and 2 on bc are lost, 2 on ce is contracted but continues outside
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l3));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;


    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));

    assert(d->getAdjList().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(d->getAdjList().front()->pl().getRouteOcc(&l2).direction == d->getAdjList().front()->getOtherNd(d));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);


    assert(tg.getNds()->size() == 4);
    // line 2 is deleted because it is lost
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, c);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    assert(tg.getNds()->size() == 4);

    // line 3 and line 2-> on bc are dead ends (because ec is contracted), line 2 on ec is also a dead end
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    assert(tg.getNds()->size() == 4);
    // all 3 lines in the triangle are lost
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);


    assert(tg.getNds()->size() == 4);

    // 3 on be and 2 on bc are lost, 2 on ce is contracted but continues outside
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l3));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);


    assert(tg.getNds()->size() == 4);
    // 1 on ec was deleted, 2 on be was deleted
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;


    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);


    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    assert(tg.getNds()->size() == 4);
    // line 2 is deleted because it is lost
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, c);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    assert(tg.getNds()->size() == 4);

    // line 3 and line 2-> on bc are dead ends (because ec is contracted), line 2 on ec is also a dead end
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);

    assert(tg.getNds()->size() == 4);
    // all 3 lines in the triangle are lost
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);


    assert(tg.getNds()->size() == 4);

    // 3 on be and 2 on bc are lost, 2 on ce is contracted but continues outside
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l3));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(e, c);


    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l1, ab, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->direction == 0);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l2, 0);

    e->pl().addConnExc(&l2, be, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    assert(tg.getNds()->size() == 4);

    // 2 on be and 2 on ce are dead ends
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->direction == 0);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, c);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, c);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, ce);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->direction == 0);

    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(d->getAdjList().front()->pl().getRouteOcc(&l2).direction == d->getAdjList().front()->getOtherNd(d));
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l3, 0);
    ce->pl().addRoute(&l2, 0);

    c->pl().addConnExc(&l2, cd, bc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->direction == 0);

    assert(a->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l2));
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().hasRoute(&l1));

    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRouteOcc(&l2).direction == 0);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRouteOcc(&l1).direction == 0);
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

    shared::transitgraph::TransitGraph tg;
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

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);
    be->pl().addRoute(&l2, 0);
    ce->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, e);

    assert(tg.getNds()->size() == 4);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->route == &l2);
    assert(tg.getEdg(a->getAdjList().front()->getOtherNd(a), d->getAdjList().front()->getOtherNd(d))->pl().getRoutes().begin()->direction == 0);
  }

  // ___________________________________________________________________________
  {
    //     1        2          1
    // a -----> b ------> c ------> d
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{100.0, 0.0}});
    auto b = tg.addNd({{200.0, 0.0}});
    auto c = tg.addNd({{300.0, 0.0}});
    auto d = tg.addNd({{400.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{100.0, 0.0}, {200.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{200.0, 0.0}, {300.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{300.0, 0.0}, {400.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(b, c);

    assert(tg.getNds()->size() == 3);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 1);

    assert(a->getAdjList().front()->pl().getRoutes().begin()->route == &l1);
    assert(d->getAdjList().front()->pl().getRoutes().begin()->route == &l1);

    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == 0);
    assert(d->getAdjList().front()->pl().getRoutes().begin()->direction == 0);
  }

  // ___________________________________________________________________________
  {
    //    1         2          1
    // a -----> b ------> c ------> d
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{100.0, 0.0}});
    auto b = tg.addNd({{200.0, 0.0}});
    auto c = tg.addNd({{300.0, 0.0}});
    auto d = tg.addNd({{400.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{100.0, 0.0}, {200.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{200.0, 0.0}, {300.0, 0.0}}});
    auto cd = tg.addEdg(c, d, {{{300.0, 0.0}, {400.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    bc->pl().addRoute(&l2, 0);
    cd->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.combineNodes(c, b);

    assert(tg.getNds()->size() == 3);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 1);

    assert(a->getAdjList().front()->pl().getRoutes().begin()->route == &l1);
    assert(d->getAdjList().front()->pl().getRoutes().begin()->route == &l1);

    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == 0);
    assert(d->getAdjList().front()->pl().getRoutes().begin()->direction == 0);
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c ---> a ---> b
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    assert(tg.getNds()->size() == 3);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ca->pl().addRoute(&l1, 0);
    ca->pl().addRoute(&l2, 0);
    ca->pl().addRoute(&l3, 0);

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    ab->pl().addRoute(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ca, ab));

    mc.combineEdges(ca, ab, a);

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

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ab, ac));

    mc.combineEdges(ab, ac, a);

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

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ac, ba));

    mc.combineEdges(ac, ba, a);

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
    // c ---> a <--- b
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    assert(tg.getNds()->size() == 3);

    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ba->pl().addRoute(&l1, 0);
    ba->pl().addRoute(&l2, 0);
    ba->pl().addRoute(&l3, 0);

    ca->pl().addRoute(&l1, 0);
    ca->pl().addRoute(&l2, 0);
    ca->pl().addRoute(&l3, 0);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ca, ba));

    mc.combineEdges(ca, ba, a);

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
    // with line directions and without connection exclusions
    // c ---> a ---> b
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    assert(tg.getNds()->size() == 3);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ca->pl().addRoute(&l1, 0);
    ca->pl().addRoute(&l2, a);
    ca->pl().addRoute(&l3, c);

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, b);
    ab->pl().addRoute(&l3, a);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ca, ab));

    mc.combineEdges(ca, ab, a);

    assert(tg.getNds()->size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    assert(b->getAdjList().front() == c->getAdjList().front());
    assert(b->getAdjList().front()->pl().getRoutes().size() == 3);

    // check geometry orientation
    assert(b->getAdjList().front()->pl().getGeom()->front().getX() == approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(b->getAdjList().front()->pl().getGeom()->back().getX() == approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == b);
      if (r.route == &l3) assert(r.direction == c);
    }
  }

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
    ab->pl().addRoute(&l2, b);
    ab->pl().addRoute(&l3, a);

    ac->pl().addRoute(&l1, 0);
    ac->pl().addRoute(&l2, a);
    ac->pl().addRoute(&l3, c);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ac, ab));

    mc.combineEdges(ab, ac, a);

    assert(tg.getNds()->size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    assert(b->getAdjList().front() == c->getAdjList().front());
    assert(b->getAdjList().front()->pl().getRoutes().size() == 3);

    // check geometry orientation
    assert(b->getAdjList().front()->pl().getGeom()->front().getX() == approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(b->getAdjList().front()->pl().getGeom()->back().getX() == approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == b);
      if (r.route == &l3) assert(r.direction == c);
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
    ba->pl().addRoute(&l2, b);
    ba->pl().addRoute(&l3, a);

    ac->pl().addRoute(&l1, 0);
    ac->pl().addRoute(&l2, a);
    ac->pl().addRoute(&l3, c);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ba, ac));

    mc.combineEdges(ac, ba, a);

    assert(tg.getNds()->size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    assert(b->getAdjList().front() == c->getAdjList().front());
    assert(b->getAdjList().front()->pl().getRoutes().size() == 3);

    // check geometry orientation
    assert(b->getAdjList().front()->pl().getGeom()->front().getX() == approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(b->getAdjList().front()->pl().getGeom()->back().getX() == approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == b);
      if (r.route == &l3) assert(r.direction == c);
    }
  }

  // ___________________________________________________________________________
  {
    // standard case without line directions and without connection exclusions
    // c ---> a <--- b
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    assert(tg.getNds()->size() == 3);

    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ba->pl().addRoute(&l1, 0);
    ba->pl().addRoute(&l2, b);
    ba->pl().addRoute(&l3, a);

    ca->pl().addRoute(&l1, 0);
    ca->pl().addRoute(&l2, a);
    ca->pl().addRoute(&l3, c);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ba, ca));

    mc.combineEdges(ca, ba, a);

    assert(tg.getNds()->size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    assert(b->getAdjList().front() == c->getAdjList().front());
    assert(b->getAdjList().front()->pl().getRoutes().size() == 3);

    // check geometry orientation
    assert(b->getAdjList().front()->pl().getGeom()->front().getX() == approx(b->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(b->getAdjList().front()->pl().getGeom()->back().getX() == approx(b->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : b->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == b);
      if (r.route == &l3) assert(r.direction == c);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c ---> a ---> b
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});

    assert(tg.getNds()->size() == 3);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 1);
    assert(c->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ca->pl().addRoute(&l1, 0);
    ca->pl().addRoute(&l2, 0);
    ca->pl().addRoute(&l3, 0);

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    ab->pl().addRoute(&l3, 0);

    a->pl().addConnExc(&l1, ca, ab);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(!mc.routeEq(ca, ab));
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c ---> a ---> b ---> d
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    assert(tg.getNds()->size() == 4);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ca->pl().addRoute(&l1, 0);
    ca->pl().addRoute(&l2, 0);
    ca->pl().addRoute(&l3, 0);

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    ab->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l3, 0);
    bd->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l3, ab, bd);
    b->pl().addConnExc(&l2, bd, ab);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ca, ab));

    mc.combineEdges(ca, ab, a);

    assert(tg.getNds()->size() == 3);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 3);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);

    // check geometry orientation
    assert(c->getAdjList().front()->pl().getGeom()->front().getX() == approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(c->getAdjList().front()->pl().getGeom()->back().getX() == approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c <--- a ---> b ---> d
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    assert(tg.getNds()->size() == 4);

    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});
    auto ab = tg.addEdg(a, b, {{{50.0, 50.0}, {100.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ac->pl().addRoute(&l1, 0);
    ac->pl().addRoute(&l2, 0);
    ac->pl().addRoute(&l3, 0);

    ab->pl().addRoute(&l1, 0);
    ab->pl().addRoute(&l2, 0);
    ab->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l3, 0);
    bd->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l3, ab, bd);
    b->pl().addConnExc(&l2, bd, ab);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ac, ab));

    mc.combineEdges(ac, ab, a);

    assert(tg.getNds()->size() == 3);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 3);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);

    // check geometry orientation
    assert(c->getAdjList().front()->pl().getGeom()->front().getX() == approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(c->getAdjList().front()->pl().getGeom()->back().getX() == approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c ---> a <--- b ---> d
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    assert(tg.getNds()->size() == 4);

    auto ca = tg.addEdg(c, a, {{{20.0, 50.0}, {50.0, 50.0}}});
    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ca->pl().addRoute(&l1, 0);
    ca->pl().addRoute(&l2, 0);
    ca->pl().addRoute(&l3, 0);

    ba->pl().addRoute(&l1, 0);
    ba->pl().addRoute(&l2, 0);
    ba->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l3, 0);
    bd->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l3, ba, bd);
    b->pl().addConnExc(&l2, bd, ba);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ba, ba));

    mc.combineEdges(ca, ba, a);

    assert(tg.getNds()->size() == 3);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 3);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);

    // check geometry orientation
    assert(c->getAdjList().front()->pl().getGeom()->front().getX() == approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(c->getAdjList().front()->pl().getGeom()->back().getX() == approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    // without line directions and with connection exclusions
    // c <--- a <--- b ---> d
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{50.0, 50.0}});
    auto b = tg.addNd({{100.0, 50.0}});
    auto c = tg.addNd({{20.0, 50.0}});
    auto d = tg.addNd({{200.0, 50.0}});

    assert(tg.getNds()->size() == 4);

    auto ac = tg.addEdg(a, c, {{{50.0, 50.0}, {20.0, 50.0}}});
    auto ba = tg.addEdg(b, a, {{{100.0, 50.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{100.0, 50.0}, {200.0, 50.0}}});

    assert(a->getAdjList().size() == 2);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ac->pl().addRoute(&l1, 0);
    ac->pl().addRoute(&l2, 0);
    ac->pl().addRoute(&l3, 0);

    ba->pl().addRoute(&l1, 0);
    ba->pl().addRoute(&l2, 0);
    ba->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l3, 0);
    bd->pl().addRoute(&l2, 0);

    b->pl().addConnExc(&l3, ba, bd);
    b->pl().addConnExc(&l2, bd, ba);

    topo::config::TopoConfig cfg;

    topo::MapConstructor mc(&cfg, &tg);

    assert(mc.routeEq(ba, ba));

    mc.combineEdges(ac, ba, a);

    assert(tg.getNds()->size() == 3);
    assert(b->getAdjList().size() == 2);
    assert(c->getAdjList().size() == 1);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 3);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 2);

    // check geometry orientation
    assert(c->getAdjList().front()->pl().getGeom()->front().getX() == approx(c->getAdjList().front()->getFrom()->pl().getGeom()->getX()));
    assert(c->getAdjList().front()->pl().getGeom()->back().getX() == approx(c->getAdjList().front()->getTo()->pl().getGeom()->getX()));

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, d
    //        c
    //     2,3^
    //   1    | 1,2,3  1,2,3
    // a ---> b ---> d ---> e

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});

    assert(tg.getNds()->size() == 5);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});

    assert(a->getAdjList().size() == 1);
    assert(b->getAdjList().size() == 3);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 2);
    assert(e->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);

    bc->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l2, 0);
    bd->pl().addRoute(&l3, 0);

    de->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l2, 0);
    de->pl().addRoute(&l3, 0);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    // should now be
    //        c
    //     2,3^
    //    1   |  1,2,3
    // a ---> d ---> e

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().size() == 1);
    assert(a->getAdjList().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 3);
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, d
    //        c
    //     2,3^
    //   1    | 1,2,3  1,2,3
    // a ---> b ---> d ---> e

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});

    assert(tg.getNds()->size() == 5);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});

    assert(a->getAdjList().size() == 1);
    assert(b->getAdjList().size() == 3);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 2);
    assert(e->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, b);

    bc->pl().addRoute(&l2, c);
    bc->pl().addRoute(&l3, c);

    bd->pl().addRoute(&l1, d);
    bd->pl().addRoute(&l2, b);
    bd->pl().addRoute(&l3, b);

    de->pl().addRoute(&l1, e);
    de->pl().addRoute(&l2, d);
    de->pl().addRoute(&l3, d);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    // should now be
    //        c
    //     2,3^
    //    1   |  1,2,3
    // a ---> d ---> e

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().size() == 1);
    assert(a->getAdjList().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 3);

    for (auto r : e->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == e);
      if (r.route == &l2) assert(r.direction == d);
      if (r.route == &l3) assert(r.direction == d);
    }

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l2) assert(r.direction == c);
      if (r.route == &l3) assert(r.direction == c);
    }

    for (auto r : a->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == d);
    }
  }

  // ___________________________________________________________________________
  {
    // node contraction of b, d
    //        c
    //     2,3^
    //   1    | 1,2,3  1,2,3
    // a ---> b ---> d ---> e

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});

    assert(tg.getNds()->size() == 5);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});

    assert(a->getAdjList().size() == 1);
    assert(b->getAdjList().size() == 3);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 2);
    assert(e->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);

    bc->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l2, 0);
    bd->pl().addRoute(&l3, 0);

    de->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l2, 0);
    de->pl().addRoute(&l3, 0);

    // no connection of line 2 from cb -> bd
    b->pl().addConnExc(&l2, bc, bd);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

    topo::config::TopoConfig cfg;
    topo::MapConstructor mc(&cfg, &tg);

    mc.combineNodes(b, d);

    // should now be
    //        c
    //     2,3^
    //    1   |  1,2,3
    // a ---> d ---> e

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().size() == 1);
    assert(a->getAdjList().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 3);

    for (auto r : e->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == 0);
      if (r.route == &l3) assert(r.direction == 0);
    }

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l2) assert(r.direction == 0);
      if (r.route == &l3) assert(r.direction == 0);
    }

    for (auto r : a->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});
    auto f = tg.addNd({{150.0, -50.0}});

    assert(tg.getNds()->size() == 6);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});
    auto df = tg.addEdg(d, f, {{{100.0, 0.0}, {150.0, -50.0}}});

    assert(a->getAdjList().size() == 1);
    assert(b->getAdjList().size() == 3);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 3);
    assert(e->getAdjList().size() == 1);
    assert(f->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);

    bc->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l2, 0);
    bd->pl().addRoute(&l3, 0);

    de->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l3, 0);

    df->pl().addRoute(&l2, 0);

    // no connection of line 3 from cb -> bd
    b->pl().addConnExc(&l3, bc, bd);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

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

    assert(tg.getNds()->size() == 5);
    assert(c->getAdjList().size() == 1);
    assert(a->getAdjList().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(f->getAdjList().front()->pl().getRoutes().size() == 1);

    for (auto r : e->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l3) assert(r.direction == 0);
    }

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l2) assert(r.direction == 0);
      if (r.route == &l3) assert(r.direction == 0);
    }

    for (auto r : a->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
    }

    for (auto r : f->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l2) assert(r.direction == 0);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{50.0, 50.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{150.0, 0.0}});
    auto f = tg.addNd({{150.0, -50.0}});

    assert(tg.getNds()->size() == 6);

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0}, {50.0, 0.0}}});
    auto bc = tg.addEdg(b, c, {{{50.0, 0.0}, {50.0, 50.0}}});
    auto bd = tg.addEdg(b, d, {{{50.0, 0.0}, {100.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {150.0, 0.0}}});
    auto df = tg.addEdg(d, f, {{{100.0, 0.0}, {150.0, -50.0}}});

    assert(a->getAdjList().size() == 1);
    assert(b->getAdjList().size() == 3);
    assert(c->getAdjList().size() == 1);
    assert(d->getAdjList().size() == 3);
    assert(e->getAdjList().size() == 1);
    assert(f->getAdjList().size() == 1);

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");

    ab->pl().addRoute(&l1, 0);

    bc->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l3, 0);

    bd->pl().addRoute(&l1, 0);
    bd->pl().addRoute(&l2, 0);
    bd->pl().addRoute(&l3, 0);

    de->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l3, 0);

    df->pl().addRoute(&l2, 0);

    // no connection of line 3 from cb -> bd
    b->pl().addConnExc(&l3, bc, bd);

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

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

    assert(tg.getNds()->size() == 5);
    assert(c->getAdjList().size() == 1);
    assert(a->getAdjList().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 2);
    assert(f->getAdjList().front()->pl().getRoutes().size() == 1);

    for (auto r : e->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l3) assert(r.direction == 0);
    }

    for (auto r : c->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l2) assert(r.direction == 0);
      if (r.route == &l3) assert(r.direction == 0);
    }

    for (auto r : a->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
    }

    for (auto r : f->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l2) assert(r.direction == 0);
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

    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 5.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {100.0, 10.0}}});
    auto ac = tg.addEdg(a, c, {{{0.0, 10.0}, {0.0, 5.0}}});
    auto cb = tg.addEdg(c, b, {{{0.0, 5.0}, {100.0, 10.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    ac->pl().addRoute(&l1, 0);
    cb->pl().addRoute(&l2, 0);

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
