// Copyright 2016
// Author: Patrick Brosi


#include <cassert>
#include <string>

#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "topo/tests/TopologicalTest2.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "topo/tests/TopoTestUtil.h"

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
void TopologicalTest2::run() {

  // ___________________________________________________________________________
  {
    //
    //                2->               2->
    //                       <----------------------------- D ^
    //          A <-------------- B                           |
    // //       |-----------------------------------------> E |
    //                             2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{250.0, 0.0}});
    auto e = tg.addNd({{500.0, 1.0}});
    auto d = tg.addNd({{500.0, 0.0}});
    // auto x = tg.addNd({{-100.0, 0.0}});
    // auto y = tg.addNd({{600.0, 0.0}});

    a->pl().addStop(transitmapper::graph::StationInfo("1", "1"));
    b->pl().addStop(transitmapper::graph::StationInfo("2", "2"));
    d->pl().addStop(transitmapper::graph::StationInfo("3", "3"));

    auto ba = tg.addEdg(b, a, {{{250.0, 5.0}, {0.0, 5.0}}});
    auto db = tg.addEdg(d, b, {{{500.0, 0.0}, {200.0, 0.0}}});


    // auto xa = tg.addEdg(x, a, {{{-100.0, 0.0}, {0.0, 0.0}}});
    // auto dy = tg.addEdg(d, y, {{{500.0, 0.0}, {600.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    // auto ae = tg.addEdg(a, e, {{{0.0, 0.0}, {500.0, 1.0}}});
    // auto ed = tg.addEdg(e, d, {{{500.0, 1.0}, {500.0, 0.0}}});
    // ae->pl().addRoute(&l2, 0);
    // ed->pl().addRoute(&l2, 0);

    // xa->pl().addRoute(&l1, 0);
    // xa->pl().addRoute(&l2, 0);

    // dy->pl().addRoute(&l1, 0);
    // dy->pl().addRoute(&l2, 0);

    ba->pl().addRoute(&l1, 0);
    ba->pl().addRoute(&l2, b);
    ba->pl().addRoute(&l1, 0);
    ba->pl().addRoute(&l2, b);
    db->pl().addRoute(&l2, d);
    db->pl().addRoute(&l1, 0);


    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, false);
    builder.removeEdgeArtifacts(&tg)                                                                                                                        ;

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    // TODO more tests

    assert(validExceptions(&tg));
    exit(1);
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //             1    |     |
    //       B ---------       -C
    //  A --------------       -D
    //           2      |     |
    //                  |-----|
    //                     2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});
    auto d = tg.addNd({{250.0, 0.0}});

    auto ad = tg.addEdg(a, d, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto bc = tg.addEdg(b, c, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    bc->pl().addRoute(&l1, 0);
    ad->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //                     1
    //                  -------
    //                  |     |
    //                  |      -C
    //  A ---B----------E      -D
    //     2     1,2    |     |
    //                  |-----|
    //                     2
    //

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);

    assert(d->getAdjList().front()->pl().hasRoute(&l2));
    assert(d->getAdjList().front()->pl().getRoutes().size() == 1);

    assert(c->getAdjList().front()->pl().hasRoute(&l1));
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //             1    |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2      |     |
    //                  |-----|
    //                     2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ac = tg.addEdg(a, c, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto bc = tg.addEdg(b, c, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ac->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //                     1
    //                  --(F)--
    //                  |     |
    //    2             |     |-C
    //  A ---B----------E     |
    //           1,2    |     |
    //                  |-(F)-|
    //                     2
    //

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != a && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == 0);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);

    assert(ec);
    assert(ef);
    assert(fc);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == 0);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == 0);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == 0);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //             1    |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2      |     |
    //                  |-----|
    //                     2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ac = tg.addEdg(a, c, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ac->pl().addRoute(&l2, 0);
    cb->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //                     1
    //                  --(F)--
    //                  |     |
    //    2             |     |-C
    //  A ---B----------E     |
    //           1,2    |     |
    //                  |-(F)-|
    //                     2
    //

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;

    for (auto nd : *tg.getNds()) if (nd->getAdjList().size() == 2 && nd->pl().getGeom()->getX() > 249) c = nd;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == 0);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);

    assert(ec);
    assert(ef);
    assert(fc);


    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == 0);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == 0);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == 0);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //             1    |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2      |     |
    //                  |-----|
    //                     2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto bc = tg.addEdg(b, c, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, 0);
    bc->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);


    //                     1
    //                  --(F)--
    //                  |     |
    //    2             |     |-C
    //  A ---B----------E     |
    //           1,2    |     |
    //                  |-(F)-|
    //                     2
    //

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == 0);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);

    assert(ec);
    assert(ef);
    assert(fc);


    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == 0);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == 0);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == 0);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }

    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //             1    |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2      |     |
    //                  |-----|
    //                     2
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, 0);
    cb->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;


    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    //                     1
    //                  --(F)--
    //                  |     |
    //    2             |     |-C
    //  A ---B----------E     |
    //           1,2    |     |
    //                  |-(F)-|
    //                     2
    //

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == 0);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);

    assert(ec);
    assert(ef);
    assert(fc);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == 0);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == 0);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == 0);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }

    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2->    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ac = tg.addEdg(a, c, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto bc = tg.addEdg(b, c, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ac->pl().addRoute(&l2, c);
    bc->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //                     1
    //                  --(F)--
    //                  |     |
    //    2->           |     |-C
    //  A ---B----------E     |
    //           1,2->  |     |
    //                  |-(F)-|
    //                     2->

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == b);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == e);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == c);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == f);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == c);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2->    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto bc = tg.addEdg(b, c, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, c);
    bc->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    //                     1
    //                  --(F)--
    //                  |     |
    //    2->           |     |-C
    //  A ---B----------E     |
    //           1,2->  |     |
    //                  |-(F)-|
    //                     2->

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == b);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == e);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == c);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == f);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == c);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2->    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, c);
    cb->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);


    //                     1
    //                  --(F)--
    //                  |     |
    //    2->           |     |-C
    //  A ---B----------E     |
    //           1,2->  |     |
    //                  |-(F)-|
    //                     2->

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == b);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == e);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == c);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == f);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == c);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           2->    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ac = tg.addEdg(a, c, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ac->pl().addRoute(&l2, c);
    cb->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    //                     1
    //                  --(F)--
    //                  |     |
    //    2->           |     |-C
    //  A ---B----------E     |
    //           1,2->  |     |
    //                  |-(F)-|
    //                     2->

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == b);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == e);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == c);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == f);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == c);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }

    assert(validExceptions(&tg));
  }

  // -----------------------------

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //          <-2     |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ac = tg.addEdg(a, c, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ac->pl().addRoute(&l2, a);
    cb->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    //                     1
    //                  --(F)--
    //                  |     |
    //    <-2           |     |-C
    //  A ---B----------E     |
    //           1,<-2  |     |
    //                  |-(F)-|
    //                    <-2

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == a);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == b);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == e);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == f);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //          <-2     |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ac = tg.addEdg(a, c, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto bc = tg.addEdg(b, c, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ac->pl().addRoute(&l2, a);
    bc->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    //                     1
    //                  --(F)--
    //                  |     |
    //    <-2           |     |-C
    //  A ---B----------E     |
    //           1,<-2  |     |
    //                  |-(F)-|
    //                    <-2

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == a);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == b);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == e);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == f);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           <-2    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto bc = tg.addEdg(b, c, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, a);
    bc->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    //                     1
    //                  --(F)--
    //                  |     |
    //    <-2           |     |-C
    //  A ---B----------E     |
    //           1,<-2  |     |
    //                  |-(F)-|
    //                    <-2

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == a);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == b);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == e);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == f);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     1
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           <-2    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, a);
    cb->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //                     1
    //                  --(F)--
    //                  |     |
    //    <-2           |     |-C
    //  A ---B----------E     |
    //           1,<-2  |     |
    //                  |-(F)-|
    //                    <-2

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == a);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l1));
    assert(be->pl().getRouteOcc(&l1).direction == 0);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == b);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l1));
      assert(ef->pl().getRouteOcc(&l1).direction == 0);

      assert(fc->pl().hasRoute(&l1));
      assert(fc->pl().getRouteOcc(&l1).direction == 0);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == e);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == f);

      assert(ec->pl().hasRoute(&l1));
      assert(ec->pl().getRouteOcc(&l1).direction == 0);
    }

    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                     <-2
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           <-2    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, a);
    cb->pl().addRoute(&l2, b);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //                    <-2
    //                  --(F)--
    //                  |     |
    //    <-2           |     |-C
    //  A ---B----------E     |
    //           <-2    |     |
    //                  |-(F)-|
    //                    <-2

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == a);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == b);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == e);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == f);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == e);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == f);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    }
    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                    2->
    //                  -------
    //                  |     |
    //       B ---------      |-C
    //  A --------------      |
    //           <-2    |     |
    //                  |-----|
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{50.0, 0.0}});
    auto c = tg.addNd({{250.0, 0.0}});

    auto ca = tg.addEdg(c, a, util::geo::PolyLine<double>({{0.0, 0.0}, {180, 0}, {180.0, -50.0}, {250.0, -50.0}, {250.0, 0.0}}));
    auto cb = tg.addEdg(c, b, util::geo::PolyLine<double>({{50.0, 0.0}, {180, 0}, {180.0, 50.0}, {250.0, 50.0}, {250.0, 0.0}}));
    ca->pl().setPolyline(ca->pl().getPolyline().getReversed());
    cb->pl().setPolyline(cb->pl().getPolyline().getReversed());

    transitmapper::graph::Route l2("2", "2", "blue");

    ca->pl().addRoute(&l2, a);
    cb->pl().addRoute(&l2, c);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    builder.removeEdgeArtifacts(&tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    //                    2->
    //                  --(F)--
    //                  |     |
    //    <-2           |     |-C
    //  A ---B----------E     |
    //           2      |     |
    //                  |-(F)-|
    //                    <-2

    b = a->getAdjList().front()->getOtherNd(a);
    TransitNode *e = 0, *f = 0;
    for (auto edg : b->getAdjList()) if (edg->getOtherNd(b) != a) e = edg->getOtherNd(b);
    for (auto edg : e->getAdjList()) if (edg->getOtherNd(e) != b && edg->getOtherNd(e) != c) f = edg->getOtherNd(e);

    assert(tg.getNds()->size() == 5);
    assert(a->getAdjList().front()->pl().hasRoute(&l2));
    assert(a->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(a->getAdjList().front()->pl().getRoutes().begin()->direction == a);

    assert(c->getAdjList().size() == 2);

    assert(e->getAdjList().size() == 3);
    assert(f->getAdjList().size() == 2);

    auto ec = tg.getEdg(e, c);
    auto ef = tg.getEdg(e, f);
    auto fc = tg.getEdg(f, c);
    auto be = tg.getEdg(b, e);

    assert(ec);
    assert(ef);
    assert(fc);
    assert(be);

    assert(be->pl().hasRoute(&l2));
    assert(be->pl().getRouteOcc(&l2).direction == 0);

    if (f->pl().getGeom()->getY() > 0) {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == f);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == c);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    } else {
      assert(ef->pl().hasRoute(&l2));
      assert(ef->pl().getRouteOcc(&l2).direction == e);

      assert(fc->pl().hasRoute(&l2));
      assert(fc->pl().getRouteOcc(&l2).direction == f);

      assert(ec->pl().hasRoute(&l2));
      assert(ec->pl().getRouteOcc(&l2).direction == e);
    }

    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //           E
    //           |
    //           |
    //           |
    //           |    1
    //           --------------------> B
    //  A --------------
    //            1    |
    //                 |
    //                 v
    //                 C
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{500.0, 0.0}});
    auto c = tg.addNd({{250.0, -200.0}});
    auto e = tg.addNd({{150.0, 200.0}});

    auto ac = tg.addEdg(a, c, util::geo::PolyLine<double>({{0.0, 0.0}, {250.0, 0.0}, {250.0, -200.0}}));
    auto eb = tg.addEdg(e, b, util::geo::PolyLine<double>({{150.0, 200.0}, {150.0, 0}, {500.0, 0.0}}));

    transitmapper::graph::Route l1("1", "1", "red");

    ac->pl().addRoute(&l1, 0);
    eb->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true, 1);
    // builder.removeEdgeArtifacts(&tg);

    // util::geo::output::GeoGraphJsonOutput gout;
    // gout.print(tg, std::cout);
    // std::cout << std::flush;

    assert(validExceptions(&tg));
  }

  // ___________________________________________________________________________
  {
    //                           2
    //                    C <------------- E
    //                    |
    //         2, 1 ->    v    1->
    //  A --------------> B --------------------> D
    //  F --------------------------------------> G
    //                     1
    //
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 0.0}});
    auto b = tg.addNd({{250.0, 0.0}});
    auto c = tg.addNd({{250.0, 0}});
    auto e = tg.addNd({{400.0, 1.0}});
    auto d = tg.addNd({{500.0, 0.0}});

    auto f = tg.addNd({{0.0, 0.0}});
    auto g = tg.addNd({{500.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 0.0},  {250.0, 0.0}}});
    auto cb = tg.addEdg(c, b, {{{250.0, 0.0},  {250.0, 0.0}}});
    auto ec = tg.addEdg(e, c, {{{400.0, 1.0},  {250.0, 0.0}}});
    auto bd = tg.addEdg(b, d, {{{250.0, 0.0},  {500.0, 0.0}}});
    auto fg = tg.addEdg(f, g, {{{0.0, 0.0},  {500.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "green");

    ab->pl().addRoute(&l1, b);
    ab->pl().addRoute(&l2, 0);
    cb->pl().addRoute(&l2, 0);
    ec->pl().addRoute(&l2, 0);
    bd->pl().addRoute(&l1, d);
    fg->pl().addRoute(&l1, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    // TODO more tests

    assert(validExceptions(&tg));
  }
}
