// Copyright 2016
// University of Freiburg - Chair of Algorithms and Datastructures
// Author: Patrick Brosi

#include <string>

#include "shared/linegraph/LineGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/tests/TopoTestUtil.h"
#include "topo/tests/TopologicalTest.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/output/GeoGraphJsonOutput.h"

#define private public
#include "topo/mapconstructor/MapConstructor.h"

using util::approx;

// _____________________________________________________________________________
void TopologicalTest::run() {
  // ___________________________________________________________________________
  {
    //     1
    // a ------> b
    // c ------> d
    //     2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 5.0}});
    auto b = tg.addNd({{50.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 5.0}, {50.0, 5.0}}});
    auto cd = tg.addEdg(c, d, {{{0.0, 0.0}, {50.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //     1, 2
    // a ------> b

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    TEST(tg.getNds()->size(), ==, 2);
    TEST((*tg.getNds()->begin())->getAdjList().front()->pl().getLines().size(),
         ==, 2);

    for (auto r :
         (*tg.getNds()->begin())->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, 0);
    }
  }
  // ___________________________________________________________________________
  {
    //      2->     1
    //     a--> b <---|
    // c -----> d <---e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto cd = tg.addEdg(c, d, {{{0.0, 0.0}, {100.0, 0.0}}});
    auto ed = tg.addEdg(e, d, {{{200.0, 0.0}, {100, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l2, b);
    cd->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    ed->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 15;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    // reconstruct the old nodes for testing
    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      }
    }

    TEST(c->getAdjList().size(), !=, 0);
    a = c->getAdjList().front()->getOtherNd(c);
    for (auto edg : a->getAdjList()) {
      if (edg->getOtherNd(a) != c) d = edg->getOtherNd(a);
    }
    for (auto edg : d->getAdjList()) {
      if (edg->getOtherNd(d) != a) e = edg->getOtherNd(d);
    }

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().lineOcc(&l2).direction, ==, c);

    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(a, d)->pl().lineOcc(&l2).direction, ==, 0);

    TEST(tg.getEdg(d, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(d, e)->pl().lineOcc(&l2).direction, ==, d);
    TEST(tg.getEdg(d, e)->pl().lineOcc(&l1).direction, ==, 0);
  }

  // ___________________________________________________________________________
  {
    //      2->     1
    //     a--> b <---|
    // c <----- d <---e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto ed = tg.addEdg(e, d, {{{200.0, 0.0}, {100, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    ed->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 15;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    // reconstruct the old nodes for testing
    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);
    for (auto edg : a->getAdjList()) {
      if (edg->getOtherNd(a) != c) d = edg->getOtherNd(a);
    }
    for (auto edg : d->getAdjList()) {
      if (edg->getOtherNd(d) != a) e = edg->getOtherNd(d);
    }

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().lineOcc(&l2).direction, ==, c);

    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(a, d)->pl().lineOcc(&l2).direction, ==, 0);

    TEST(tg.getEdg(d, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(d, e)->pl().lineOcc(&l2).direction, ==, d);
    TEST(tg.getEdg(d, e)->pl().lineOcc(&l1).direction, ==, 0);
  }

  // ___________________________________________________________________________
  {
    //      2->     1
    //     a--> b <---|
    // c <----- d --->e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    de->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 15;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);
    for (auto edg : a->getAdjList()) {
      if (edg->getOtherNd(a) != c) d = edg->getOtherNd(a);
    }
    for (auto edg : d->getAdjList()) {
      if (edg->getOtherNd(d) != a) e = edg->getOtherNd(d);
    }

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().lineOcc(&l2).direction, ==, c);

    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);

    TEST(tg.getEdg(a, d)->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(a, d)->pl().lineOcc(&l2).direction, ==, 0);

    TEST(tg.getEdg(d, e)->pl().getLines().size(), ==, 2);
    TEST(tg.getEdg(d, e)->pl().lineOcc(&l2).direction, ==, d);
    TEST(tg.getEdg(d, e)->pl().lineOcc(&l1).direction, ==, 0);
  }

  // ___________________________________________________________________________
  {
    //             1
    //          b <---|
    // c -----> d <---e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto cd = tg.addEdg(d, c, {{{0.0, 0.0}, {100.0, 0.0}}});
    auto ed = tg.addEdg(e, d, {{{200, 0.0}, {100.0, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    cd->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    ed->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    // reconstruct the old nodes for testing
    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      }
    }

    cd = c->getAdjList().front();
    for (auto* nd : *tg.getNds()) {
      if (nd->getDeg() == 1 && nd != c) e = nd;
      if (nd->getDeg() == 2) d = nd;
    }
    ed = e->getAdjList().front();
    TEST(cd->pl().hasLine(&l2) && cd->pl().lineOcc(&l2).direction, ==, c);
    TEST(ed->pl().hasLine(&l2) && ed->pl().lineOcc(&l2).direction, !=, 0);
    TEST(ed->pl().lineOcc(&l2).direction, !=, e);
  }

  // ___________________________________________________________________________
  {
    //             1
    //          b <---|
    // c <----- d --->e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    dc->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    de->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //    <-2          1,<-2
    // c ---------> d < --- e

    // reconstruct the old nodes for testing
    c = 0;
    e = 0;
    d = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      } else {
        d = nd;
      }
    }

    TEST(tg.getNds()->size(), ==, 3);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);

    for (auto ro : e->getAdjList().front()->pl().getLines()) {
      if (ro.line == &l1) TEST(ro.direction, ==, 0);
      if (ro.line == &l2) TEST(ro.direction, ==, d);
    }
  }
  // ___________________________________________________________________________
  {
    //             1
    //          b <---|
    // c <----- d <---e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto ed = tg.addEdg(e, d, {{{200, 0.0}, {100.0, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    dc->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    ed->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();

    //    <-2          1,<-2
    // c ---------> d < --- e

    // reconstruct the old nodes for testing
    c = 0;
    e = 0;
    d = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      } else {
        d = nd;
      }
    }

    TEST(tg.getNds()->size(), ==, 3);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);

    for (auto ro : e->getAdjList().front()->pl().getLines()) {
      if (ro.line == &l1) TEST(ro.direction, ==, 0);
      if (ro.line == &l2) TEST(ro.direction, ==, d);
    }
  }

  // ___________________________________________________________________________
  {
    //      2->     1
    //     a<-- b <---|
    // c <------ d --->e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{96.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ba = tg.addEdg(b, a, {{{95.0, 10.0}, {30.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {95, 10.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ba->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    de->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 20;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //    <-2   2        1,<-2
    // c ----a-----> d < --- e

    // reconstruct the old nodes for testing
    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      }
    }

    assert(c && e);
    a = c->getAdjList().front()->getOtherNd(c);
    d = e->getAdjList().front()->getOtherNd(e);

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);

    for (auto edg : a->getAdjList()) {
      if (edg->getOtherNd(a) == c) continue;
      TEST(edg->pl().getLines().size(), ==, 1);
      TEST(edg->pl().getLines().begin()->direction, ==, 0);
    }

    for (auto ro : e->getAdjList().front()->pl().getLines()) {
      if (ro.line == &l1) TEST(ro.direction, ==, 0);
      if (ro.line == &l2) TEST(ro.direction, ==, d);
    }
  }

  // ___________________________________________________________________________
  {
    //     1
    // a ------> b
    // c <------ d
    //     2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{50.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {50.0, 10.0}}});
    auto cd = tg.addEdg(d, c, {{{50.0, 0.0}, {0.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 15;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();

    //     1, 2
    // a ------> b

    TEST(tg.getNds()->size(), ==, 2);
    TEST((*tg.getNds()->begin())->getAdjList().front()->pl().getLines().size(),
         ==, 2);

    for (auto r :
         (*tg.getNds()->begin())->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    //     1
    // a ------> b
    // c <------ d
    //     2->
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{50.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {50.0, 10.0}}});
    auto cd = tg.addEdg(d, c, {{{50.0, 0.0}, {0.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 15;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //     1, 2
    // a ------> b

    TEST(tg.getNds()->size(), ==, 2);
    TEST((*tg.getNds()->begin())->getAdjList().front()->pl().getLines().size(),
         ==, 2);

    for (auto r :
         (*tg.getNds()->begin())->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) {
        TEST(r.direction, !=, 0);
        TEST(r.direction->pl().getGeom()->getX(), >, 25);
      }
    }
  }

  // ___________________________________________________________________________
  {
    //     1
    // a ----------> b
    // c <---------- d
    //      2->
    // e <---------- f
    //    <-3,2->
    // g <---------- h
    //      3
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{0.0, 5.0}});
    auto b = tg.addNd({{50.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto e = tg.addNd({{0.0, 10.0}});
    auto f = tg.addNd({{50.0, 10.0}});
    auto g = tg.addNd({{0.0, 15.0}});
    auto h = tg.addNd({{50.0, 15.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 5.0}, {50.0, 5.0}}});
    auto cd = tg.addEdg(d, c, {{{50.0, 0.0}, {0.0, 0.0}}});

    auto fe = tg.addEdg(f, e, {{{50.0, 10.0}, {0.0, 10.0}}});
    auto hg = tg.addEdg(h, g, {{{50.0, 15.0}, {0.0, 15.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");
    shared::linegraph::Line l3("3", "3", "green");
    shared::linegraph::Line l4("4", "4", "black");

    ab->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, d);
    fe->pl().addLine(&l3, e);
    fe->pl().addLine(&l2, f);
    hg->pl().addLine(&l3, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 20;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //  1, 2->, 3
    // a ------> b

    TEST(tg.getNds()->size(), ==, 2);
    TEST((*tg.getNds()->begin())->getAdjList().front()->pl().getLines().size(),
         ==, 3);

    for (auto r :
         (*tg.getNds()->begin())->getAdjList().front()->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) {
        TEST(r.direction, !=, 0);
        TEST(r.direction->pl().getGeom()->getX(), >, 25);
      }
      if (r.line == &l3) TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    //      1
    //    a->b
    // c ------> d
    //     2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{70.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 5.0}, {70.0, 5.0}}});
    auto cd = tg.addEdg(c, d, {{{0.0, 0.0}, {100.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    cd->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //    2   1,2    2
    // c ----a--->b----> d
    c = 0;
    d = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (!c)
          c = nd;
        else
          d = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);
    b = d->getAdjList().front()->getOtherNd(d);

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(a, b), !=, 0);
    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 2);

    for (auto r : tg.getEdg(a, b)->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    //      1
    //     a->b
    // c <----- d
    //     2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{70.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 5.0}, {70.0, 5.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    dc->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //    2   1,2    2
    // c ----a--->b----> d

    c = 0;
    d = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (!c)
          c = nd;
        else
          d = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);
    b = d->getAdjList().front()->getOtherNd(d);

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(d->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(tg.getEdg(a, b), !=, 0);
    TEST(tg.getEdg(a, b)->pl().getLines().size(), ==, 2);

    for (auto r : tg.getEdg(a, b)->pl().getLines()) {
      if (r.line == &l1) TEST(r.direction, ==, 0);
      if (r.line == &l2) TEST(r.direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    //      1
    //     a--> b
    // c <----- d
    //     2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{100.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 5.0}, {100.0, 5.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l1, 0);
    dc->pl().addLine(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //    2     1,2
    // c ----a-----> d
    //
    c = 0;
    d = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          d = nd;
        else
          c = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);

    TEST(tg.getNds()->size(), ==, 3);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      TEST(e->pl().getLines().size(), ==, 2);

      for (auto r : e->pl().getLines()) {
        if (r.line == &l1) TEST(r.direction, ==, 0);
        if (r.line == &l2) TEST(r.direction, ==, 0);
      }
    }
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a--> b
    // c <----- d
    //     <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{100.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 5.0}, {100.0, 5.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});

    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //    <-2   2
    // c ----a-----> d
    //
    c = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1 && nd->pl().getGeom()->getX() < 50) {
        c = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);

    TEST(tg.getNds()->size(), ==, 3);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      TEST(e->pl().getLines().size(), ==, 1);
      TEST(e->pl().getLines().begin()->direction, ==, 0);
    }
  }

  // ___________________________________________________________________________
  {
    //      2->     1
    //     a--> b <---|
    // c <----- d --->e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{100.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 5.0}, {100.0, 5.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 5.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);
    eb->pl().addLine(&l1, 0);
    de->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();

    //    <-2   2        1,<-2
    // c ----a-----> d < --- e

    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l1))
          e = nd;
        else
          c = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 2);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      TEST(e->pl().getLines().size(), ==, 1);
      TEST(e->pl().getLines().begin()->direction, ==, 0);
    }

    for (auto ro : e->getAdjList().front()->pl().getLines()) {
      if (ro.line == &l1) TEST(ro.direction, ==, 0);
      if (ro.line == &l2)
        TEST(ro.direction, ==, e->getAdjList().front()->getOtherNd(e));
    }
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a--> b
    // c <----- d --->e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{100.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 5.0}, {100.0, 5.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);
    de->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();

    //    <-2   2        <-2
    // c ----a-----> d < --- e

    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l2)) {
          if (nd->getAdjList().front()->pl().lineOcc(&l2).direction != nd)
            e = nd;
          else
            c = nd;
        }
      }
    }

    assert(c && e);

    a = c->getAdjList().front()->getOtherNd(c);
    d = e->getAdjList().front()->getOtherNd(e);

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 1);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      TEST(e->pl().getLines().size(), ==, 1);
      TEST(e->pl().getLines().begin()->direction, ==, 0);
    }

    TEST(e->getAdjList().front()->pl().getLines().begin()->direction, ==, d);
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a--> b
    // c <----- d <---e
    //     <-2    <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{100.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 5.0}, {100.0, 5.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto ed = tg.addEdg(e, d, {{{200.0, 0.0}, {100, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ab->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);
    ed->pl().addLine(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();

    //    <-2   2        <-2
    // c ----a-----> d < --- e

    c = 0;
    e = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l2)) {
          if (nd->getAdjList().front()->pl().lineOcc(&l2).direction != nd)
            e = nd;
          else
            c = nd;
        }
      }
    }

    assert(c && e);

    a = c->getAdjList().front()->getOtherNd(c);
    d = e->getAdjList().front()->getOtherNd(e);

    TEST(tg.getNds()->size(), ==, 4);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);
    TEST(e->getAdjList().front()->pl().getLines().size(), ==, 1);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      TEST(e->pl().getLines().size(), ==, 1);
      TEST(e->pl().getLines().begin()->direction, ==, 0);
    }

    TEST(e->getAdjList().front()->pl().getLines().begin()->direction, ==, d);
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a<-- b
    // c <----- d
    //     <-2
    shared::linegraph::LineGraph tg;
    auto a = tg.addNd({{30.0, 5.0}});
    auto b = tg.addNd({{100.0, 5.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ba = tg.addEdg(b, a, {{{100.0, 5.0}, {30.0, 5.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});

    shared::linegraph::Line l1("1", "1", "red");
    shared::linegraph::Line l2("2", "2", "blue");

    ba->pl().addLine(&l2, b);
    dc->pl().addLine(&l2, c);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 10;

    topo::MapConstructor mc(&cfg, &tg);
    mc.collapseShrdSegs();
    mc.removeEdgeArtifacts();

    //    <-2   2
    // c ----a-----> d

    c = 0;
    d = 0;
    for (auto nd : *tg.getNds()) {
      if (nd->getDeg() == 1) {
        if (nd->getAdjList().front()->pl().hasLine(&l2)) {
          if (nd->getAdjList().front()->pl().lineOcc(&l2).direction == nd)
            c = nd;
          else
            d = nd;
        }
      }
    }

    assert(c && d);

    a = c->getAdjList().front()->getOtherNd(c);

    TEST(tg.getNds()->size(), ==, 3);
    TEST(c->getAdjList().front()->pl().getLines().size(), ==, 1);
    TEST(c->getAdjList().front()->pl().getLines().begin()->direction, ==, c);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      TEST(e->pl().getLines().size(), ==, 1);
      TEST(e->pl().getLines().begin()->direction, ==, 0);
    }
  }
}
