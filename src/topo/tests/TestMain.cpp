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
    //     1
    // a ------> b
    // c ------> d
    //     2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{50.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {50.0, 10.0}}});
    auto cd = tg.addEdg(c, d, {{{0.0, 0.0}, {50.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //     1, 2
    // a ------> b

    assert(tg.getNds()->size() == 2);
    assert((*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes().size() == 2);

    for (auto r : (*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    //     1
    // a ------> b
    // c <------ d
    //     2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{50.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {50.0, 10.0}}});
    auto cd = tg.addEdg(d, c, {{{50.0, 0.0}, {0.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //     1, 2
    // a ------> b

    assert(tg.getNds()->size() == 2);
    assert((*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes().size() == 2);

    for (auto r : (*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    //     1
    // a ------> b
    // c <------ d
    //     2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{50.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {50.0, 10.0}}});
    auto cd = tg.addEdg(d, c, {{{50.0, 0.0}, {0.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //     1, 2
    // a ------> b

    assert(tg.getNds()->size() == 2);
    assert((*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes().size() == 2);

    for (auto r : (*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) {
        assert(r.direction != 0);
        assert(r.direction->pl().getGeom()->getX() == approx(50));
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
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{0.0, 10.0}});
    auto b = tg.addNd({{50.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{50.0, 0.0}});

    auto e = tg.addNd({{0.0, 20.0}});
    auto f = tg.addNd({{50.0, 20.0}});
    auto g = tg.addNd({{0.0, 30.0}});
    auto h = tg.addNd({{50.0, 30.0}});

    auto ab = tg.addEdg(a, b, {{{0.0, 10.0}, {50.0, 10.0}}});
    auto cd = tg.addEdg(d, c, {{{50.0, 0.0}, {0.0, 0.0}}});

    auto fe = tg.addEdg(f, e, {{{50.0, 20.0}, {0.0, 20.0}}});
    auto hg = tg.addEdg(h, g, {{{50.0, 30.0}, {0.0, 30.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");
    transitmapper::graph::Route l3("3", "3", "green");
    transitmapper::graph::Route l4("4", "4", "black");

    ab->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, d);
    fe->pl().addRoute(&l3, e);
    fe->pl().addRoute(&l2, f);
    hg->pl().addRoute(&l3, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //  1, 2->, 3
    // a ------> b

    assert(tg.getNds()->size() == 2);
    assert((*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes().size() == 3);

    for (auto r : (*tg.getNds()->begin())->getAdjList().front()->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) {
        assert(r.direction != 0);
        assert(r.direction->pl().getGeom()->getX() == approx(50));
      }
      if (r.route == &l3) {
        assert(r.direction == 0);
      }
    }
  }

  // ___________________________________________________________________________
  {
    //      1
    //    a->b
    // c ------> d
    //     2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{70.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {70.0, 10.0}}});
    auto cd = tg.addEdg(c, d, {{{0.0, 0.0}, {100.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    cd->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);


    //    2   1,2    2
    // c ----a--->b----> d

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a, b) != 0);
    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 2);

    for (auto r : tg.getEdg(a, b)->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    //      1
    //     a->b
    // c <----- d
    //     2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{70.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {70.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    dc->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);


    //    2   1,2    2
    // c ----a--->b----> d

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(d->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(tg.getEdg(a, b) != 0);
    assert(tg.getEdg(a, b)->pl().getRoutes().size() == 2);

    for (auto r : tg.getEdg(a, b)->pl().getRoutes()) {
      if (r.route == &l1) assert(r.direction == 0);
      if (r.route == &l2) assert(r.direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    //      1
    //     a--> b
    // c <----- d
    //     2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l1, 0);
    dc->pl().addRoute(&l2, 0);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //    2     1,2
    // c ----a-----> d

    assert(tg.getNds()->size() == 3);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      assert(e->pl().getRoutes().size() == 2);

      for (auto r : e->pl().getRoutes()) {
        if (r.route == &l1) assert(r.direction == 0);
        if (r.route == &l2) assert(r.direction == 0);
      }
    }
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a--> b
    // c <----- d
    //     <-2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});

    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l2, b);
    dc->pl().addRoute(&l2, c);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //    <-2   2
    // c ----a-----> d

    assert(tg.getNds()->size() == 3);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().begin()->direction == c);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      assert(e->pl().getRoutes().size() == 1);
      assert(e->pl().getRoutes().begin()->direction == 0);
    }
  }

  // ___________________________________________________________________________
  {
    //      2->     1
    //     a--> b <---|
    // c <----- d --->e
    //     <-2    <-2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l2, b);
    dc->pl().addRoute(&l2, c);
    eb->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l2, d);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);

    //    <-2   2        1,<-2
    // c ----a-----> d < --- e

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().begin()->direction == c);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 2);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      assert(e->pl().getRoutes().size() == 1);
      assert(e->pl().getRoutes().begin()->direction == 0);
    }

    for (auto ro : e->getAdjList().front()->pl().getRoutes()) {
      if (ro.route == &l1) assert(ro.direction == 0);
      if (ro.route == &l2) assert(ro.direction->pl().getGeom()->getX() == approx(100));
    }
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a--> b
    // c <----- d --->e
    //     <-2    <-2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l2, b);
    dc->pl().addRoute(&l2, c);
    de->pl().addRoute(&l2, d);

    d->pl().addConnExc(&l2, de, dc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);



    //    <-2   2        <-2
    // c ----a-----> d < --- e

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().begin()->direction == c);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 1);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      assert(e->pl().getRoutes().size() == 1);
      assert(e->pl().getRoutes().begin()->direction == 0);
    }

    assert(e->getAdjList().front()->pl().getRoutes().begin()->direction->pl().getGeom()->getX() == approx(100));
    assert(e->getAdjList().front()->pl().getRoutes().begin()->direction->pl().connOccurs(&l2, e->getAdjList().front(), c->getAdjList().front()));
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a--> b
    // c <----- d <---e
    //     <-2    <-2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto ed = tg.addEdg(e, d, {{{200.0, 0.0}, {100, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l2, b);
    dc->pl().addRoute(&l2, c);
    ed->pl().addRoute(&l2, d);

    d->pl().addConnExc(&l2, ed, dc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    builder.removeEdgeArtifacts(&tg);


    //    <-2   2        <-2
    // c ----a-----> d < --- e

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().begin()->direction == c);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 1);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      assert(e->pl().getRoutes().size() == 1);
      assert(e->pl().getRoutes().begin()->direction == 0);
    }

    assert(e->getAdjList().front()->pl().getRoutes().begin()->direction->pl().getGeom()->getX() == approx(100));
    assert(e->getAdjList().front()->pl().getRoutes().begin()->direction->pl().connOccurs(&l2, e->getAdjList().front(), c->getAdjList().front()));
  }

  // ___________________________________________________________________________
  {
    //      2->
    //     a<-- b
    // c <----- d <---e
    //     <-2    <-2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ba = tg.addEdg(a, b, {{{100.0, 10.0}, {30.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto ed = tg.addEdg(e, d, {{{200.0, 0.0}, {100, 0.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ba->pl().addRoute(&l2, b);
    dc->pl().addRoute(&l2, c);
    ed->pl().addRoute(&l2, d);

    d->pl().addConnExc(&l2, ed, dc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    std::cerr << "START" << std::endl;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    // builder.removeEdgeArtifacts(&tg);

    util::geo::output::GeoGraphJsonOutput gout;
    gout.print(tg, std::cout);

    //    <-2   2        <-2
    // c ----a-----> d < --- e

    // update the node variables, may have changed
    for (auto nd : *tg.getNds()) {
      if (nd->getAdjList().size() == 1) {
        if (nd->getAdjList().front()->pl().getRoutes().begin()->direction == nd) {
          c = nd;
        }
        else e = nd;
      }
    }

    a = c->getAdjList().front()->getOtherNd(c);

    assert(a->pl().getGeom()->getX() == approx(30));

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().begin()->direction == c);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 1);

    for (auto e : a->getAdjList()) {
      std::cerr << e->getFrom() << " -> " << e->getTo() << std::endl;
      if (e->getOtherNd(a) == c) continue;
      assert(e->pl().getRoutes().size() == 1);
      assert(e->pl().getRoutes().begin()->direction == 0);
    }

    assert(e->getAdjList().front()->pl().getRoutes().begin()->direction->pl().getGeom()->getX() == approx(100));
    assert(e->getAdjList().front()->pl().getRoutes().begin()->direction->pl().connOccurs(&l2, e->getAdjList().front(), c->getAdjList().front()));
  }

  // ___________________________________________________________________________
  {
    //      2->     1
    //     a--> b <---|
    // c <----- d --->e
    //     <-2    <-2
    shared::transitgraph::TransitGraph tg;
    auto a = tg.addNd({{30.0, 10.0}});
    auto b = tg.addNd({{100.0, 10.0}});
    auto c = tg.addNd({{0.0, 0.0}});
    auto d = tg.addNd({{100.0, 0.0}});
    auto e = tg.addNd({{200.0, 0.0}});

    auto ab = tg.addEdg(a, b, {{{30.0, 10.0}, {100.0, 10.0}}});
    auto dc = tg.addEdg(d, c, {{{100.0, 0.0}, {0.0, 0.0}}});
    auto de = tg.addEdg(d, e, {{{100.0, 0.0}, {200, 0.0}}});
    auto eb = tg.addEdg(e, b, {{{200.0, 0.0}, {100, 10.0}}});

    transitmapper::graph::Route l1("1", "1", "red");
    transitmapper::graph::Route l2("2", "2", "blue");

    ab->pl().addRoute(&l2, b);
    dc->pl().addRoute(&l2, c);
    eb->pl().addRoute(&l1, 0);
    de->pl().addRoute(&l2, d);

    d->pl().addConnExc(&l2, de, dc);

    topo::config::TopoConfig cfg;
    cfg.maxAggrDistance = 50;

    topo::Builder builder(&cfg);
    builder.createTopologicalNodes(&tg, true);
    // builder.removeEdgeArtifacts(&tg);


    //    <-2   2        1,<-2
    // c ----a-----> d < --- e

    assert(tg.getNds()->size() == 4);
    assert(c->getAdjList().front()->pl().getRoutes().size() == 1);
    assert(c->getAdjList().front()->pl().getRoutes().begin()->direction == c);
    assert(e->getAdjList().front()->pl().getRoutes().size() == 2);

    for (auto e : a->getAdjList()) {
      if (e->getOtherNd(a) == c) continue;
      assert(e->pl().getRoutes().size() == 1);
      assert(e->pl().getRoutes().begin()->direction == 0);
    }

    for (auto ro : e->getAdjList().front()->pl().getRoutes()) {
      if (ro.route == &l1) assert(ro.direction == 0);
      if (ro.route == &l2) assert(ro.direction->pl().getGeom()->getX() == approx(100));
    }
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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ca, ab));

    builder.combineEdges(ca, ab, a, &tg);

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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ab, ac));

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

    assert(builder.routeEq(ac, ba));

    builder.combineEdges(ac, ba, a, &tg);

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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ca, ba));

    builder.combineEdges(ca, ba, a, &tg);

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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ca, ab));

    builder.combineEdges(ca, ab, a, &tg);

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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ac, ab));

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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ba, ac));

    builder.combineEdges(ac, ba, a, &tg);

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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ba, ca));

    builder.combineEdges(ca, ba, a, &tg);

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

    topo::Builder builder(&cfg);

    assert(!builder.routeEq(ca, ab));
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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ca, ab));

    builder.combineEdges(ca, ab, a, &tg);

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

    assert(!b->pl().connOccurs(&l3, c->getAdjList().front(), bd));
    assert(!b->pl().connOccurs(&l2, c->getAdjList().front(), bd));

    assert(!b->pl().connOccurs(&l3, bd, c->getAdjList().front()));
    assert(!b->pl().connOccurs(&l2, bd, c->getAdjList().front()));
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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ac, ab));

    builder.combineEdges(ac, ab, a, &tg);

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

    assert(!b->pl().connOccurs(&l3, c->getAdjList().front(), bd));
    assert(!b->pl().connOccurs(&l2, c->getAdjList().front(), bd));

    assert(!b->pl().connOccurs(&l3, bd, c->getAdjList().front()));
    assert(!b->pl().connOccurs(&l2, bd, c->getAdjList().front()));
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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ba, ba));

    builder.combineEdges(ca, ba, a, &tg);

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

    assert(!b->pl().connOccurs(&l3, c->getAdjList().front(), bd));
    assert(!b->pl().connOccurs(&l2, c->getAdjList().front(), bd));

    assert(!b->pl().connOccurs(&l3, bd, c->getAdjList().front()));
    assert(!b->pl().connOccurs(&l2, bd, c->getAdjList().front()));
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

    topo::Builder builder(&cfg);

    assert(builder.routeEq(ba, ba));

    builder.combineEdges(ac, ba, a, &tg);

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

    assert(!b->pl().connOccurs(&l2, c->getAdjList().front(), bd));
    assert(!b->pl().connOccurs(&l3, c->getAdjList().front(), bd));

    assert(!b->pl().connOccurs(&l3, bd, c->getAdjList().front()));
    assert(!b->pl().connOccurs(&l2, bd, c->getAdjList().front()));
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
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

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
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

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

    assert(!b->pl().connOccurs(&l2, bc, bd));

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

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

    //
    auto contred = e->getAdjList().front()->getFrom();

    assert(contred->pl().getConnExc().size() == 1);
    assert(contred->pl().getConnExc().begin()->first == &l2);

    assert(!contred->pl().connOccurs(&l2, c->getAdjList().front(), e->getAdjList().front()));
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

    assert(!b->pl().connOccurs(&l3, bc, bd));

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(b, d, &tg);

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


    auto contred = e->getAdjList().front()->getFrom();

    assert(contred->pl().getConnExc().size() == 1);
    assert(contred->pl().getConnExc().begin()->first == &l3);

    assert(!contred->pl().connOccurs(&l3, c->getAdjList().front(), e->getAdjList().front()));
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

    assert(!b->pl().connOccurs(&l3, bc, bd));

    assert(c->getAdjList().front()->pl().getRoutes().size() == 2);

    topo::config::TopoConfig cfg;
    topo::Builder builder(&cfg);

    builder.combineNodes(d, b, &tg);

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
    auto contred = e->getAdjList().front()->getFrom();

    assert(contred->pl().getConnExc().size() == 1);
    assert(contred->pl().getConnExc().begin()->first == &l3);

    assert(!contred->pl().connOccurs(&l3, c->getAdjList().front(), e->getAdjList().front()));
  }
}
