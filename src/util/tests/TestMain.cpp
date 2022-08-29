// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include <clocale>
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/tests/QuadTreeTest.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/graph/Algorithm.h"
#include "util/graph/Dijkstra.h"
#include "util/graph/BiDijkstra.h"
#include "util/graph/DirGraph.h"
#include "util/graph/EDijkstra.h"
#include "util/graph/UndirGraph.h"
#include "util/json/Writer.h"

using namespace util;
using namespace util::geo;
using namespace util::graph;

using util::approx;

// _____________________________________________________________________________
int main(int argc, char** argv) {
	UNUSED(argc);
	UNUSED(argv);

  QuadTreeTest quadTreeTest;
  quadTreeTest.run();

  // ___________________________________________________________________________
  {
    TEST(geo::frechetDist(Line<double>{{0, 0}, {10, 10}}, Line<double>{{0, 0},
                                  {10, 10}}, 1) == approx(0));

    TEST(geo::frechetDist(Line<double>{{0, 0}, {10, 10}}, Line<double>{{0, 0},
                                  {10, 10}}, 0.1) == approx(0));

    TEST(geo::frechetDist(Line<double>{{0, 10}, {10, 10}}, Line<double>{{0, 0},
                                  {10, 0}}, 0.1) == approx(10));

    TEST(geo::frechetDist(Line<double>{{0, 10}, {10, 10}}, Line<double>{{0, 0},
                                  {0, 0}}, 0.1)
        == approx(util::geo::dist(Point<double>{0, 0}, Point<double>{10, 10})));

    TEST(geo::frechetDist(Line<double>{{0, 10}, {0, 10}}, Line<double>{{0, 0},
                                  {0, 5}}, 0.1) == approx(10));

    TEST(geo::frechetDist(Line<double>{{0, 10}, {0, 10}}, Line<double>{{0, 5},
                                  {0, 5}}, 0.1) == approx(5));
  }

  // ___________________________________________________________________________
  {
    TEST(util::btsSimi("", ""), ==, approx(1));
    TEST(util::btsSimi("Hallo", "Test"), ==, approx(0));
    TEST(util::btsSimi("Test", "Hallo"), ==, approx(0));
    TEST(util::btsSimi("Test", "Test"), ==, approx(1));
    TEST(util::btsSimi("Milner Road / Wandlee Road", "Wandlee Road"), ==, approx(1));
    TEST(util::btsSimi("bla blubb blob", "blubb blib"), ==, approx(0.9));
    TEST(util::btsSimi("St Pancras International", "London St Pancras"), ==, approx(0.588235));
    TEST(util::btsSimi("Reiterstrasse", "Reiterstrasse Freiburg im Breisgau"), ==, approx(1));
    TEST(util::btsSimi("Reiterstrasse", "Reiter Freiburg im Breisgau"), ==, approx(.53333333));
    TEST(util::btsSimi("AA", "Reiterstrasse, Freiburg im Breisgau"), ==, approx(0));
    TEST(util::btsSimi("blibb blabbel bla blubb blob", "blubb blib blabb"), ==, approx(0.875));
    TEST(util::btsSimi("blibb blabbel bla blubb blobo", "blubb blib blabb blabo"), ==, approx(0.84));
    TEST(util::btsSimi("blubb blib blabb", "blibb blabbel bla blubb blob"), ==, approx(0.875));
    TEST(util::btsSimi("blubbb blib blabb blobo", "blibb blabbel bla blubb blobo"), ==, approx(0.84));
    TEST(util::btsSimi("Reiter Freiburg im Breisgau", "Reiter Frei burg im Brei sgau"), ==, approx(0.931034));
    // fallback to jaccard
    TEST(util::btsSimi("Freiburg im Breisgau, Germany, Main Railway Station", "Main Railway Station Freiburg im Breisgau, Germany"), ==, approx(1));

  }

  // ___________________________________________________________________________
  {
    std::string test = u8"Zuerich, Hauptbahnhof (Nord)";
    auto tokens = util::tokenize(test);

    TEST(tokens.size(), ==, 3);

    TEST(util::jaccardSimi("Zuerich Hauptbahnhof Nord", "Zuerich, Hauptbahnhof (Nord)"), ==, approx(1));
    TEST(util::jaccardSimi("Zuerich Hauptbahnhof", "Zuerich, Hauptbahnhof ()"), ==, approx(1));
    TEST(util::jaccardSimi("Zuerich Hauptbahnhof", "Zuerich, Hauptbahnhof (Nord)"), ==, approx(2./3.));
  }

  // ___________________________________________________________________________
  {
    TEST(util::atof("45.534215"), ==, approx(45.534215));
    TEST(util::atof("5.534"), ==, approx(5.534));
    TEST(util::atof("534"), ==, approx(534));
    TEST(util::atof("-534"), ==, approx(-534));
    TEST(util::atof("-45.534215"), ==, approx(-45.534215));
    TEST(util::atof("-45.534215", 2), ==, approx(-45.53));

    // TODO: more test cases
  }

  // ___________________________________________________________________________
  {
    std::stringstream ss;
    util::json::Writer wr(&ss, 2, false);

    util::json::Val a("bla");
    util::json::Val b(1);
    util::json::Val c(1.0);
    util::json::Val d("a");
    util::json::Val e({"a", "b", "c"});

    util::json::Val f({1, json::Array{2, 3, 4}, 3});

    ss.str("");
    wr = util::json::Writer(&ss, 2, false);
    util::json::Val i({1, json::Array{2, json::Null(), 4}, true});
    wr.val(i);
    wr.closeAll();
    TEST(ss.str(), ==, "[1,[2,null,4],true]");

    ss.str("");
    wr = util::json::Writer(&ss, 2, false);
    i = util::json::Val({1, json::Array{2, json::Null(), 4}, false});
    wr.val(i);
    wr.closeAll();
    TEST(ss.str(), ==, "[1,[2,null,4],false]");

    ss.str("");
    wr = util::json::Writer(&ss, 2, false);
    i = util::json::Val({1, json::Array{2, json::Null(), 4}, false});
    wr.val(i);
    wr.closeAll();
    TEST(ss.str(), ==, "[1,[2,null,4],false]");

    ss.str("");
    wr = util::json::Writer(&ss, 2, false);
    i = util::json::Val({1, json::Array{2.13, "", 4}, 0});
    wr.val(i);
    wr.closeAll();
    TEST(ss.str(), ==, "[1,[2.13,\"\",4],0]");

    ss.str("");
    wr = util::json::Writer(&ss, 2, false);
    i = util::json::Val(
        {1, json::Array{2.13, json::Dict{{"a", 1}, {"B", 2.123}}, 4}, 0});
    wr.val(i);
    wr.closeAll();
    TEST((ss.str() == "[1,[2.13,{\"a\":1,\"B\":2.12},4],0]" ||
            ss.str() == "[1,[2.13,{\"B\":2.12,\"a\":1},4],0]"));
  }

  // ___________________________________________________________________________
  {
    DirGraph<int, int> g;

    DirNode<int, int>* a = new DirNode<int, int>(0);
    DirNode<int, int>* b = new DirNode<int, int>(0);
    g.addNd(a);
    TEST(g.getNds().size(), ==, (size_t)1);
    g.addNd(b);
    TEST(g.getNds().size(), ==, (size_t)2);

    g.addEdg(a, b);
    TEST(a->getDeg(), ==, (size_t)1);
    TEST(b->getDeg(), ==, (size_t)0);

    auto c = g.addNd();

    g.addEdg(a, c);
    g.addEdg(c, b);
    TEST(a->getDeg(), ==, (size_t)2);
    TEST(b->getDeg(), ==, (size_t)0);
    TEST(c->getDeg(), ==, (size_t)1);

    g.delEdg(a, c);

    TEST(a->getDeg(), ==, (size_t)1);
    TEST(b->getDeg(), ==, (size_t)0);
    TEST(c->getDeg(), ==, (size_t)1);

    g.addEdg(a, a);
    TEST(a->getDeg(), ==, (size_t)2);

    g.delEdg(a, a);
    TEST(a->getDeg(), ==, (size_t)1);

    g.delEdg(a, a);
    TEST(a->getDeg(), ==, (size_t)1);

    // TODO: more test cases
  }

  // ___________________________________________________________________________
  {
    UndirGraph<int, int> g;

    UndirNode<int, int>* a = new UndirNode<int, int>(0);
    UndirNode<int, int>* b = new UndirNode<int, int>(0);
    g.addNd(a);
    TEST(g.getNds().size(), ==, (size_t)1);
    g.addNd(b);
    TEST(g.getNds().size(), ==, (size_t)2);

    g.addEdg(a, b);
    TEST(a->getDeg(), ==, (size_t)1);
    TEST(b->getDeg(), ==, (size_t)1);

    auto c = g.addNd();

    g.addEdg(a, c);
    g.addEdg(c, b);
    TEST(a->getDeg(), ==, (size_t)2);
    TEST(b->getDeg(), ==, (size_t)2);
    TEST(c->getDeg(), ==, (size_t)2);

    g.delEdg(a, c);

    TEST(a->getDeg(), ==, (size_t)1);
    TEST(b->getDeg(), ==, (size_t)2);
    TEST(c->getDeg(), ==, (size_t)1);

    g.delNd(b);

    TEST(a->getDeg(), ==, (size_t)0);
    TEST(c->getDeg(), ==, (size_t)0);

    g.addEdg(a, a);
    TEST(a->getDeg(), ==, (size_t)1);

    g.delEdg(a, a);
    TEST(a->getDeg(), ==, (size_t)0);

    // TODO: more test cases
  }

  // ___________________________________________________________________________
  {
    Grid<int, Line, double> g(
        .5, .5, Box<double>(Point<double>(0, 0), Point<double>(3, 3)));

    Line<double> l;
    l.push_back(Point<double>(0, 0));
    l.push_back(Point<double>(1.5, 2));

    Line<double> l2;
    l2.push_back(Point<double>(2.5, 1));
    l2.push_back(Point<double>(2.5, 2));

    g.add(l, 1);
    g.add(l2, 2);

    std::set<int> ret;

    Box<double> req(Point<double>(.5, 1), Point<double>(1, 1.5));
    g.get(req, &ret);
    TEST(ret.size(), ==, (size_t)1);

    ret.clear();
    g.getNeighbors(1, 0, &ret);
    TEST(ret.size(), ==, (size_t)1);

    ret.clear();
    g.getNeighbors(1, 0.55, &ret);
    TEST(ret.size(), ==, (size_t)2);

    // TODO: more test cases
  }

  // ___________________________________________________________________________
  {
    Line<double> a;
    a.push_back(Point<double>(1, 1));
    a.push_back(Point<double>(10, 1));

    auto dense = util::geo::densify(a, 1);

    TEST(dense.size(), ==, (size_t)10);

    for (int i = 0; i < 10; i++) {
      TEST(dense[i].getX(), ==, approx(i + 1.0));
    }

    dense = util::geo::simplify(dense, 0.1);
    TEST(dense.size(), ==, (size_t)2);

    Line<double> b;
    b.push_back(Point<double>(1, 1));
    b.push_back(Point<double>(5, 7));
    b.push_back(Point<double>(10, 3));

    dense = util::geo::densify(b, 1);

    dense = util::geo::simplify(dense, 0.1);
    TEST(dense.size(), ==, (size_t)3);
  }

  // ___________________________________________________________________________
  {
    Line<double> a;
    a.push_back(Point<double>(1, 1));
    a.push_back(Point<double>(2, 1));
    a.push_back(Point<double>(3, 1));
    a.push_back(Point<double>(3, 2));
    a.push_back(Point<double>(4, 2));
    a.push_back(Point<double>(4, 1));
    a.push_back(Point<double>(5, 1));
    a.push_back(Point<double>(6, 1));

    Line<double> b;
    b.push_back(Point<double>(1, 1));
    b.push_back(Point<double>(2, 1));
    b.push_back(Point<double>(3, 1));
    b.push_back(Point<double>(4, 1));
    b.push_back(Point<double>(5, 1));
    b.push_back(Point<double>(6, 1));

    double fd = util::geo::accFrechetDistC(a, b, 0.1);
    TEST(fd, ==, approx(2));
  }

  // ___________________________________________________________________________
  {
    Line<double> e;
    e.push_back(Point<double>(1, 1));
    e.push_back(Point<double>(1, 2));

    Line<double> f;
    f.push_back(Point<double>(1, 1));
    f.push_back(Point<double>(1, 2));

    double fd = util::geo::frechetDist(e, f, 0.1);

    TEST(fd, ==, approx(0));

    Line<double> a;
    a.push_back(Point<double>(1, 1));
    a.push_back(Point<double>(2, 1));
    a.push_back(Point<double>(3, 2));
    a.push_back(Point<double>(4, 2));
    a.push_back(Point<double>(5, 1));
    a.push_back(Point<double>(6, 1));

    Line<double> b;
    b.push_back(Point<double>(1, 1));
    b.push_back(Point<double>(2, 1));
    b.push_back(Point<double>(3, 1));
    b.push_back(Point<double>(4, 1));
    b.push_back(Point<double>(5, 1));
    b.push_back(Point<double>(6, 1));

    auto adense = util::geo::densify(a, 0.1);
    auto bdense = util::geo::densify(b, 0.1);

    fd = util::geo::frechetDist(a, b, 0.1);

    TEST(fd, ==, approx(1));

    Line<double> c;
    c.push_back(Point<double>(1, 1));
    c.push_back(Point<double>(2, 1));

    Line<double> d;
    d.push_back(Point<double>(3, 1));
    d.push_back(Point<double>(4, 1));

    fd = util::geo::frechetDist(c, d, 0.1);

    TEST(fd, ==, approx(2));

    Line<double> g;
    g.push_back(Point<double>(1, 1));
    g.push_back(Point<double>(10, 1));

    Line<double> h;
    h.push_back(Point<double>(1, 1));
    h.push_back(Point<double>(3, 2));
    h.push_back(Point<double>(3, 1));
    h.push_back(Point<double>(10, 1));

    fd = util::geo::frechetDist(g, h, 0.1);

    TEST(fd, ==, approx(1));
  }

  // ___________________________________________________________________________
  {
    Line<double> a;
    a.push_back(Point<double>(1, 1));
    a.push_back(Point<double>(1, 2));

    Line<double> b;
    b.push_back(Point<double>(1, 2));
    b.push_back(Point<double>(2, 2));

    Line<double> c;
    c.push_back(Point<double>(2, 2));
    c.push_back(Point<double>(2, 1));

    Line<double> d;
    d.push_back(Point<double>(2, 1));
    d.push_back(Point<double>(1, 1));

    Box<double> box(Point<double>(2, 3), Point<double>(5, 4));
    MultiLine<double> ml;
    ml.push_back(a);
    ml.push_back(b);
    ml.push_back(c);
    ml.push_back(d);

    TEST(parallelity(box, ml), ==, approx(1));
    ml = rotate(ml, 45);
    TEST(parallelity(box, ml), ==, approx(0));
    ml = rotate(ml, 45);
    TEST(parallelity(box, ml), ==, approx(1));
    ml = rotate(ml, 45);
    TEST(parallelity(box, ml), ==, approx(0));
    ml = rotate(ml, 45);
    TEST(parallelity(box, ml), ==, approx(1));
  }

  // ___________________________________________________________________________
  {
    TEST("zürich", ==, util::urlDecode("z%C3%BCrich"));
    TEST("!@$%^*()", ==, util::urlDecode("!%40%24%25%5E*()"));
    TEST("Løkken", ==, util::urlDecode("L%C3%B8kken"));
    TEST("á é", ==, util::urlDecode("%C3%A1%20%C3%A9"));
    TEST("á é", ==, util::urlDecode("%C3%A1+%C3%A9"));
  }

  // ___________________________________________________________________________
  {
    TEST("Hello\\\\Goodbye!" == util::jsonStringEscape("Hello\\Goodbye!"));
    TEST("\\\"Hello\\\"" == util::jsonStringEscape("\"Hello\""));
  }

  // ___________________________________________________________________________
  {
    TEST(util::split("hello,again", ',').size(), ==, (size_t)2);
    TEST(util::split("hello,,again", ',').size(), ==, (size_t)3);
    TEST(util::split("hello", ',').size(), ==, (size_t)1);
    TEST(util::split("", ',').size(), ==, (size_t)0);
  }

  // ___________________________________________________________________________
  {
    TEST(util::editDist("hello", "mello"), ==, (size_t)1);
    TEST(util::editDist("mello", "hello"), ==, (size_t)1);
    TEST(util::editDist("abcde", "abfde"), ==, (size_t)1);
    TEST(util::editDist("abcd", "abcde"), ==, (size_t)1);
    TEST(util::editDist("xabcd", "abcde"), ==, (size_t)2);
    TEST(util::editDist("abcd", "abcdes"), ==, (size_t)2);
    TEST(util::editDist("hello", "hello"), ==, (size_t)0);
  }

  // ___________________________________________________________________________
  {
    TEST(util::prefixEditDist("hello", "hello", 0), ==, (size_t)0);
    TEST(util::prefixEditDist("hello", "hello", 100), ==, (size_t)0);
    TEST(util::prefixEditDist("hello", "hello"), ==, (size_t)0);
    TEST(util::prefixEditDist("hel", "hello"), ==, (size_t)0);
    TEST(util::prefixEditDist("hel", "hello", 0), ==, (size_t)0);
    TEST(util::prefixEditDist("hel", "hello", 1), ==, (size_t)0);
    TEST(util::prefixEditDist("hel", "hello", 2), ==, (size_t)0);
    TEST(util::prefixEditDist("hal", "hello", 2), ==, (size_t)1);
    TEST(util::prefixEditDist("hal", "hello", 1), ==, (size_t)1);
    TEST(util::prefixEditDist("hal", "hello", 0), >, (size_t)0);
    TEST(util::prefixEditDist("fel", "hello", 0), >, (size_t)0);
    TEST(util::prefixEditDist("fel", "hello", 1), ==, (size_t)1);
    TEST(util::prefixEditDist("fel", "hello", 2), ==, (size_t)1);
    TEST(util::prefixEditDist("fal", "hello", 2), ==, (size_t)2);
    TEST(util::prefixEditDist("fal", "hello", 1), >, (size_t)1);
    TEST(util::prefixEditDist("fal", "hello", 0), >, (size_t)0);
    TEST(util::prefixEditDist("far", "hello", 0), >, (size_t)0);
    TEST(util::prefixEditDist("far", "hello", 1), >, (size_t)1);
    TEST(util::prefixEditDist("far", "hello", 2), >, (size_t)2);
    TEST(util::prefixEditDist("far", "hello", 3), ==, (size_t)3);
    TEST(util::prefixEditDist("far", "hello", 4), ==, (size_t)3);
    TEST(util::prefixEditDist("far", "hello"), ==, (size_t)3);
    TEST(util::prefixEditDist("hefar", "hello"), ==, (size_t)3);
    TEST(util::prefixEditDist("hefaree", "hello"), ==, (size_t)5);
    TEST(util::prefixEditDist("helloo", "hello"), ==, (size_t)1);
    TEST(util::prefixEditDist("helloo", "hello", 0), >, (size_t)0);
    TEST(util::prefixEditDist("helloo", "hello", 1), ==, (size_t)1);
    TEST(util::prefixEditDist("helloo", "hello", 2), ==, (size_t)1);
    TEST(util::prefixEditDist("", "hello", 2), ==, (size_t)0);
    TEST(util::prefixEditDist("e", "hello", 2), ==, (size_t)1);
    TEST(util::prefixEditDist("el", "hello", 2), ==, (size_t)1);
    TEST(util::prefixEditDist("ello", "hello", 2), ==, (size_t)1);
    TEST(util::prefixEditDist("hell", "hello", 2), ==, (size_t)0);
    TEST(util::prefixEditDist("hell", "", 2), >, (size_t)2);
    TEST(util::prefixEditDist("hell", ""), ==, (size_t)4);
  }

  // ___________________________________________________________________________
  {
    TEST(util::toString(34) == "34");
    TEST(util::toString("34") == "34");
  }

  // ___________________________________________________________________________
  {
    std::string a("lorem ipsum ipsum lorem");

    TEST(util::replace(a, "ips", "aa"));
    TEST(a, ==, "lorem aaum ipsum lorem");

    TEST(!util::replace(a, "blablabla", ""));
    TEST(a, ==, "lorem aaum ipsum lorem");

    TEST(util::replace(a, "m", ""));
    TEST(a, ==, "lore aaum ipsum lorem");

    TEST(!util::replace(a, "", ""));
    TEST(a, ==, "lore aaum ipsum lorem");

    std::string b("lorem ipsum ipsum lorem");
    TEST(util::replaceAll(b, "ips", "aa"));
    TEST(b, ==, "lorem aaum aaum lorem");

    TEST(util::replaceAll(b, "m", ""));
    TEST(b, ==, "lore aau aau lore");

    TEST(util::replaceAll(b, "a", "aa"));
    TEST(b, ==, "lore aaaau aaaau lore");

    TEST(util::replaceAll(b, "e", "e"));
    TEST(b, ==, "lore aaaau aaaau lore");

    TEST(util::replaceAll(b, "e", "ee"));
    TEST(b, ==, "loree aaaau aaaau loree");

    TEST(!util::replaceAll(b, "", "ee"));
    TEST(b, ==, "loree aaaau aaaau loree");
  }

  // ___________________________________________________________________________
  {
    UndirGraph<std::string, int> g;

    auto a = g.addNd("A");
    auto b = g.addNd("B");
    auto c = g.addNd("C");
    auto d = g.addNd("D");
    auto e = g.addNd("E");

    g.addEdg(a, c, 1);
    g.addEdg(a, b, 5);
    g.addEdg(d, c, 1);
    g.addEdg(d, b, 3);
    g.addEdg(e, d, 1);
    g.addEdg(e, b, 1);

    auto comps = util::graph::Algorithm::connectedComponents(g);

    TEST(comps.size(), ==, static_cast<size_t>(1));
    TEST(comps[0].count(a));
    TEST(comps[0].count(b));
    TEST(comps[0].count(c));
    TEST(comps[0].count(d));
    TEST(comps[0].count(e));

    auto f = g.addNd("F");
    comps = util::graph::Algorithm::connectedComponents(g);
    TEST(comps.size(), ==, static_cast<size_t>(2));

    auto gn = g.addNd("G");
    comps = util::graph::Algorithm::connectedComponents(g);
    TEST(comps.size(), ==, static_cast<size_t>(3));

    g.addEdg(f, gn, 1);
    comps = util::graph::Algorithm::connectedComponents(g);
    TEST(comps.size(), ==, static_cast<size_t>(2));

    g.addEdg(f, a, 1);
    comps = util::graph::Algorithm::connectedComponents(g);
    TEST(comps.size(), ==, static_cast<size_t>(1));
  }

  // ___________________________________________________________________________
  {
    DirGraph<std::string, int> g;

    auto a = g.addNd("A");
    auto b = g.addNd("B");
    auto c = g.addNd("C");
    auto d = g.addNd("D");
    auto e = g.addNd("E");

    auto eAC = g.addEdg(a, c, 1);
    auto eAB = g.addEdg(a, b, 5);
    auto eDC = g.addEdg(d, c, 1);
    auto eDB = g.addEdg(d, b, 3);
    auto eED = g.addEdg(e, d, 1);
    auto eEB = g.addEdg(e, b, 1);

    UNUSED(eAC);
    UNUSED(eDC);
    UNUSED(eDB);
    UNUSED(eED);
    UNUSED(eEB);

    struct CostFunc : public EDijkstra::CostFunc<std::string, int, int> {
      int operator()(const Edge<std::string, int>* from,
                     const Node<std::string, int>* n,
                     const Edge<std::string, int>* to) const {
        UNUSED(from);

        // dont count cost of start edge
        if (n) return to->pl();
        return 0;
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    auto cost = EDijkstra::shortestPath(eAB, cFunc);

    for (auto u : cost) {
      int single = EDijkstra::shortestPath(eAB, u.first, cFunc);
      TEST(single, ==, u.second);
    }

    // all to 1
    auto eBC = g.addEdg(b, c, 10);

    auto costb = EDijkstra::shortestPathRev(eBC, cFunc);
    for (auto u : costb) {
      int single = EDijkstra::shortestPath(u.first, eBC, cFunc);
      TEST(single, ==, u.second);
    }
  }

  // ___________________________________________________________________________
  {
    UndirGraph<std::string, int> g;

    auto a = g.addNd("A");
    auto b = g.addNd("B");
    auto c = g.addNd("C");
    auto d = g.addNd("D");
    auto e = g.addNd("E");

    auto eAC = g.addEdg(a, c, 1);
    auto eAB = g.addEdg(a, b, 5);
    auto eDC = g.addEdg(d, c, 1);
    auto eDB = g.addEdg(d, b, 3);
    auto eED = g.addEdg(e, d, 1);
    auto eEB = g.addEdg(e, b, 1);

    UNUSED(eAC);
    UNUSED(eDC);
    UNUSED(eDB);
    UNUSED(eED);
    UNUSED(eEB);

    struct CostFunc : public EDijkstra::CostFunc<std::string, int, int> {
      int operator()(const Edge<std::string, int>* from,
                     const Node<std::string, int>* n,
                     const Edge<std::string, int>* to) const {
        UNUSED(from);

        // dont count cost of start edge
        if (n) return to->pl();
        return 0;
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    EDijkstra::NList<std::string, int> res;
    EDijkstra::EList<std::string, int> resE;
    int cost = EDijkstra::shortestPath(eAB, d, cFunc, &resE, &res);

    TEST(cost, ==, 2);

    TEST(resE.size(), ==, (size_t)3);
    TEST(res.size(), ==, (size_t)3);
    TEST((*(res.rbegin()))->pl(), ==, "A");
    TEST((*(++res.rbegin()))->pl(), ==, "C");
    TEST((*(++++res.rbegin()))->pl(), ==, "D");

    TEST((*(resE.rbegin())), ==, eAB);
    TEST((*(++resE.rbegin())), ==, eAC);
    TEST((*(++++resE.rbegin())), ==, eDC);

    cost = EDijkstra::shortestPath(eAB, b, cFunc, &resE, &res);
    TEST(cost, ==, 0);
  }

  // ___________________________________________________________________________
  {
    UndirGraph<std::string, int> g;

    auto a = g.addNd("A");
    auto b = g.addNd("B");
    auto c = g.addNd("C");
    auto d = g.addNd("D");
    auto e = g.addNd("E");

    auto eAC = g.addEdg(a, c, 1);
    auto eAB = g.addEdg(a, b, 5);
    auto eDC = g.addEdg(d, c, 1);
    auto eDB = g.addEdg(d, b, 3);
    auto eED = g.addEdg(e, d, 1);
    auto eEB = g.addEdg(e, b, 1);

    UNUSED(eAC);
    UNUSED(eDC);
    UNUSED(eDB);
    UNUSED(eED);
    UNUSED(eEB);

    struct CostFunc : public EDijkstra::CostFunc<std::string, int, int> {
      int operator()(const Edge<std::string, int>* from,
                     const Node<std::string, int>* n,
                     const Edge<std::string, int>* to) const {
        UNUSED(from);

        // dont count cost of start edge
        if (n) return to->pl();
        return 0;
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    std::set<Node<std::string, int>*> tos;
    tos.insert(d);
    tos.insert(b);
    tos.insert(b);

    EDijkstra::NList<std::string, int> res;
    EDijkstra::EList<std::string, int> resE;
    int cost = EDijkstra::shortestPath(eAB, tos, cFunc, &resE, &res);
    TEST(cost, ==, 0);
  }

  // ___________________________________________________________________________
  {
    UndirGraph<std::string, int> g;

    auto a = g.addNd("A");
    auto b = g.addNd("B");
    auto c = g.addNd("C");
    auto d = g.addNd("D");
    auto e = g.addNd("E");

    g.addEdg(a, c, 1);
    auto eAB = g.addEdg(a, b, 5);
    auto eDC = g.addEdg(d, c, 1);
    g.addEdg(d, b, 3);
    auto eED = g.addEdg(e, d, 1);
    g.addEdg(e, b, 1);

    struct CostFunc : public EDijkstra::CostFunc<std::string, int, int> {
      int operator()(const Edge<std::string, int>* from,
                     const Node<std::string, int>* n,
                     const Edge<std::string, int>* to) const {
        UNUSED(from);

        // dont count cost of start edge
        if (n) return to->pl();
        return 0;
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    std::set<Edge<std::string, int>*> tos;
    tos.insert(eDC);
    tos.insert(eED);

    std::unordered_map<Edge<std::string, int>*,
                       EDijkstra::EList<std::string, int>*>
        resE;
    resE[eDC] = new EDijkstra::EList<std::string, int>();
    resE[eED] = new EDijkstra::EList<std::string, int>();
    std::unordered_map<Edge<std::string, int>*,
                       EDijkstra::NList<std::string, int>*>
        res;
    res[eDC] = new EDijkstra::NList<std::string, int>();
    res[eED] = new EDijkstra::NList<std::string, int>();
    auto hFunc = ZeroHeurFunc<std::string, int, int>();
    std::unordered_map<Edge<std::string, int>*, int> cost =
        EDijkstra::shortestPath(eAB, tos, cFunc, hFunc, resE, res);

    TEST(cost[eDC], ==, 2);
    TEST(cost[eED], ==, 2);

    TEST(resE[eDC]->size(), ==, (size_t)3);
    TEST(res[eED]->size(), ==, (size_t)3);

    TEST(resE[eDC]->size(), ==, (size_t)3);
    TEST(res[eED]->size(), ==, (size_t)3);
  }

  // ___________________________________________________________________________
  {
    UndirGraph<std::string, int> g;

    auto a = g.addNd("A");
    auto b = g.addNd("B");
    auto c = g.addNd("C");
    auto d = g.addNd("D");
    auto e = g.addNd("E");

    g.addEdg(a, c, 1);
    g.addEdg(a, b, 5);
    g.addEdg(d, c, 1);
    g.addEdg(d, b, 3);
    g.addEdg(e, d, 1);
    g.addEdg(e, b, 1);

    struct CostFunc : public EDijkstra::CostFunc<std::string, int, int> {
      int operator()(const Edge<std::string, int>* fr,
                     const Node<std::string, int>* n,
                     const Edge<std::string, int>* to) const {
        UNUSED(fr);
        UNUSED(n);
        return to->pl();
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    EDijkstra::NList<std::string, int> res;
    EDijkstra::EList<std::string, int> resE;
    int cost = EDijkstra::shortestPath(a, b, cFunc, &resE, &res);

    TEST(res.size(), ==, (size_t)5);
    TEST((*(res.rbegin()))->pl(), ==, "A");
    TEST((*(++res.rbegin()))->pl(), ==, "C");
    TEST((*(++++res.rbegin()))->pl(), ==, "D");
    TEST((*(++++++res.rbegin()))->pl(), ==, "E");
    TEST((*(++++++++res.rbegin()))->pl(), ==, "B");
    TEST(cost, ==, 4);
    TEST((*(resE.rbegin()))->getFrom()->pl(), ==, "A");
    TEST((*(++resE.rbegin()))->getFrom()->pl(), ==, "D");
    TEST((*(++++resE.rbegin()))->getFrom()->pl(), ==, "E");
    TEST((*(++++++resE.rbegin()))->getTo()->pl(), ==, "B");

    TEST(resE.size(), ==, (size_t)4);

    cost = EDijkstra::shortestPath(d, b, cFunc, &res);
    TEST(cost, ==, 2);

    cost = EDijkstra::shortestPath(b, d, cFunc, &res);
    TEST(cost, ==, 2);

    cost = EDijkstra::shortestPath(e, b, cFunc, &res);
    TEST(cost, ==, 1);

    cost = EDijkstra::shortestPath(b, e, cFunc, &res);
    TEST(cost, ==, 1);

    cost = EDijkstra::shortestPath(b, a, cFunc, &res);
    TEST(cost, ==, 4);

    cost = EDijkstra::shortestPath(c, a, cFunc, &res);
    TEST(cost, ==, 1);

    cost = EDijkstra::shortestPath(a, c, cFunc, &res);
    TEST(cost, ==, 1);

    cost = EDijkstra::shortestPath(a, d, cFunc, &res);
    TEST(cost, ==, 2);
  }

  // ___________________________________________________________________________
  {
    DirGraph<int, int> g;

    DirNode<int, int>* a = new DirNode<int, int>(1);
    DirNode<int, int>* b = new DirNode<int, int>(4);
    g.addNd(a);
    g.addNd(b);

    auto c = g.addNd(2);
    auto d = g.addNd(3);
    auto x = g.addNd();

    g.addEdg(a, d, 4);
    g.addEdg(a, c, 1);
    g.addEdg(c, b, 1);
    g.addEdg(b, d, 1);

    struct CostFunc : public EDijkstra::CostFunc<int, int, int> {
      int operator()(const Edge<int, int>* fr, const Node<int, int>* n,
                     const Edge<int, int>* to) const {
        UNUSED(fr);
        UNUSED(n);
        return to->pl();
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    EDijkstra::NList<int, int> res;
    int cost = EDijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 3);

    g.addEdg(c, d, 3);
    cost = EDijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 3);

    g.addEdg(a, b, 1);
    g.addEdg(x, a, 1);
    cost = EDijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 2);
  }

  // ___________________________________________________________________________
  {
    DirGraph<int, int> g;

    auto source = g.addNd();
    auto target = g.addNd();
    auto a = g.addNd();
    auto b = g.addNd();

    g.addEdg(source, a, 4);
    g.addEdg(source, b, 5);
    g.addEdg(a, target, 3);
    g.addEdg(b, target, 1);

    struct CostFunc : public BiDijkstra::CostFunc<int, int, int> {
      int operator()(const Node<int, int>* fr, const Edge<int, int>* e,
                     const Node<int, int>* to) const {
        UNUSED(fr);
        UNUSED(to);
        return e->pl();
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    BiDijkstra::NList<int, int> res;
    int cost = BiDijkstra::shortestPath(source, target, cFunc, &res);

    TEST(cost, ==, 6);
  }

  // ___________________________________________________________________________
  {
    DirGraph<int, int> g;

    DirNode<int, int>* a = new DirNode<int, int>(1);
    DirNode<int, int>* b = new DirNode<int, int>(0);
    g.addNd(a);
    g.addNd(b);

    auto c = g.addNd();
    auto d = g.addNd(4);
    auto x = g.addNd();

    g.addEdg(a, d, 4);
    g.addEdg(a, c, 1);
    g.addEdg(c, b, 1);
    g.addEdg(b, d, 1);

    struct CostFunc : public BiDijkstra::CostFunc<int, int, int> {
      int operator()(const Node<int, int>* fr, const Edge<int, int>* e,
                     const Node<int, int>* to) const {
        UNUSED(fr);
        UNUSED(to);
        return e->pl();
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    BiDijkstra::NList<int, int> res;
    int cost = BiDijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 3);
    TEST(res.size(), ==, (size_t)4);

    g.addEdg(c, d, 3);
    cost = BiDijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 3);

    g.addEdg(a, b, 1);
    g.addEdg(x, a, 1);
    cost = BiDijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 2);

    // const std::set<Node<int, int>*> to{b, c, d, x};
    // std::unordered_map<Node<int, int>*, BiDijkstra::EList<int, int>*> resEdges;
    // std::unordered_map<Node<int, int>*, BiDijkstra::NList<int, int>*> resNodes;

    // for (auto n : to) {
      // resEdges[n] = new BiDijkstra::EList<int, int>();
      // resNodes[n] = new BiDijkstra::NList<int, int>();
    // }

    // auto costs = BiDijkstra::shortestPath(a, to, cFunc, resEdges, resNodes);

    // TEST(costs[b], ==, 1);
    // TEST(costs[c], ==, 1);
    // TEST(costs[d], ==, 2);
    // TEST(costs[x], ==, 999);
  }

  // ___________________________________________________________________________
  {
    DirGraph<int, int> g;

    DirNode<int, int>* a = new DirNode<int, int>(1);
    DirNode<int, int>* b = new DirNode<int, int>(0);
    g.addNd(a);
    g.addNd(b);

    auto c = g.addNd();
    auto d = g.addNd(4);
    auto x = g.addNd();

    g.addEdg(a, d, 4);
    g.addEdg(a, c, 1);
    g.addEdg(c, b, 1);
    g.addEdg(b, d, 1);

    struct CostFunc : public Dijkstra::CostFunc<int, int, int> {
      int operator()(const Node<int, int>* fr, const Edge<int, int>* e,
                     const Node<int, int>* to) const {
        UNUSED(fr);
        UNUSED(to);
        return e->pl();
      };
      int inf() const { return 999; };
    };

    CostFunc cFunc;

    Dijkstra::NList<int, int> res;
    int cost = Dijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 3);
    TEST(res.size(), ==, (size_t)4);

    g.addEdg(c, d, 3);
    cost = Dijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 3);

    g.addEdg(a, b, 1);
    g.addEdg(x, a, 1);
    cost = Dijkstra::shortestPath(a, d, cFunc, &res);

    TEST(cost, ==, 2);

    const std::set<Node<int, int>*> to{b, c, d, x};
    std::unordered_map<Node<int, int>*, Dijkstra::EList<int, int>*> resEdges;
    std::unordered_map<Node<int, int>*, Dijkstra::NList<int, int>*> resNodes;

    for (auto n : to) {
      resEdges[n] = new Dijkstra::EList<int, int>();
      resNodes[n] = new Dijkstra::NList<int, int>();
    }

    auto costs = Dijkstra::shortestPath(a, to, cFunc, resEdges, resNodes);

    TEST(costs[b], ==, 1);
    TEST(costs[c], ==, 1);
    TEST(costs[d], ==, 2);
    TEST(costs[x], ==, 999);
  }

  // ___________________________________________________________________________
  {{util::Nullable<std::string> nullable;
  TEST(nullable.isNull());
}

{
  util::Nullable<std::string> nullable(0);
  TEST(nullable.isNull());
}

{
  std::string str = "aa";
  util::Nullable<std::string> nullable(&str);
  TEST(!nullable.isNull());

  TEST(nullable == "aa");
  TEST(!(nullable == "aaa"));
  TEST(!(nullable != "aa"));
  TEST(nullable == "aa");

  TEST(nullable.get(), ==, "aa");
  TEST(std::string(nullable), ==, "aa");
}

{
  int a = 23;
  util::Nullable<int> nullable(a);
  util::Nullable<int> nullable2(24);
  TEST(!nullable.isNull());

  TEST(nullable, ==, 23);
  TEST(nullable, >=, 23);
  TEST(nullable, <=, 23);
  TEST(nullable, <, 24);
  TEST(nullable, <, 24);
  TEST(!(nullable < 22));
  TEST(nullable, !=, nullable2);
  TEST(nullable, <, nullable2);
  TEST(nullable2, >, nullable);

  util::Nullable<int> nullable3(nullable);
  TEST(nullable == nullable3);

  nullable3 = nullable2;
  TEST(nullable2 == nullable3);
  TEST(nullable3 == 24);
  TEST(nullable2 == 24);
  TEST(nullable2 == nullable2.get());
  TEST(int(nullable2) == nullable2.get());
  TEST(!nullable3.isNull());
  TEST(!nullable2.isNull());

  util::Nullable<int> voidnull;
  TEST(voidnull.isNull());
}
}

// ___________________________________________________________________________
{
  auto p = pointFromWKT<double>("POINT(10 50)");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("POINT( 10 50)");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("POINT (10 50 30)");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("POINT (10     50 30)");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("POINT(10 50 30)");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("POINT (10    50) ");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("MPOINT(10 50 30)");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("MPOINT(10 50)");
  TEST(p.getX(), ==, approx(10));
  TEST(p.getY(), ==, approx(50));

  p = pointFromWKT<double>("POINT(10.05 50.05)");
  TEST(p.getX(), ==, approx(10.05));
  TEST(p.getY(), ==, approx(50.05));

  auto wktl = lineFromWKT<double>("LINESTRING(0 0, 1 1,2 3, 0 1)");
  TEST(wktl.size(), ==, (size_t)4);
  TEST(wktl[0].getX(), ==, approx(0));
  TEST(wktl[0].getY(), ==, approx(0));
  TEST(wktl[1].getX(), ==, approx(1));
  TEST(wktl[1].getY(), ==, approx(1));
  TEST(wktl[2].getX(), ==, approx(2));
  TEST(wktl[2].getY(), ==, approx(3));
  TEST(wktl[3].getX(), ==, approx(0));
  TEST(wktl[3].getY(), ==, approx(1));

  wktl = lineFromWKT<double>("MLINESTRING(0 0, 1 1,2 3, 0 1)");
  TEST(wktl.size(), ==, (size_t)4);
  TEST(wktl[0].getX(), ==, approx(0));
  TEST(wktl[0].getY(), ==, approx(0));
  TEST(wktl[1].getX(), ==, approx(1));
  TEST(wktl[1].getY(), ==, approx(1));
  TEST(wktl[2].getX(), ==, approx(2));
  TEST(wktl[2].getY(), ==, approx(3));
  TEST(wktl[3].getX(), ==, approx(0));
  TEST(wktl[3].getY(), ==, approx(1));

  wktl = lineFromWKT<double>("MLINESTRING (0 0, 1  1,2   3, 0 1 )");
  TEST(wktl.size(), ==, (size_t)4);
  TEST(wktl[0].getX(), ==, approx(0));
  TEST(wktl[0].getY(), ==, approx(0));
  TEST(wktl[1].getX(), ==, approx(1));
  TEST(wktl[1].getY(), ==, approx(1));
  TEST(wktl[2].getX(), ==, approx(2));
  TEST(wktl[2].getY(), ==, approx(3));
  TEST(wktl[3].getX(), ==, approx(0));
  TEST(wktl[3].getY(), ==, approx(1));
}
// ___________________________________________________________________________
{
  geo::Point<double> a(1, 2);
  geo::Point<double> b(2, 3);
  geo::Point<double> c(4, 5);
  TEST(a.getX(), ==, approx(1));
  TEST(a.getY(), ==, approx(2));

  a.setX(3);
  TEST(a.getX(), ==, approx(3));
  TEST(a.getY(), ==, approx(2));

  a.setY(4);
  TEST(a.getX(), ==, approx(3));
  TEST(a.getY(), ==, approx(4));

  auto d = a + b;
  TEST(d.getX(), ==, approx(5));
  TEST(d.getY(), ==, approx(7));

  a.setX(1);
  a.setY(2);

  TEST(geo::dist(a, a), ==, approx(0));
  TEST(geo::dist(a, b), ==, approx(sqrt(2)));

  d = d + d;

  geo::Box<double> box(a, c);
  TEST(geo::contains(a, box));
  TEST(geo::contains(b, box));
  TEST(geo::contains(c, box));
  TEST(!geo::contains(d, box));

  geo::Line<double> line{a, b, c};

  TEST(geo::contains(line, box));
  line.push_back(d);
  TEST(!geo::contains(line, box));

  geo::LineSegment<double> ls{a, b};
  TEST(geo::contains(a, ls));
  TEST(geo::contains(b, ls));
  TEST(!geo::contains(c, ls));
  TEST(geo::contains(a + geo::Point<double>(.5, .5), ls));
  TEST(!geo::contains(a + geo::Point<double>(1.5, 1.5), ls));

  geo::LineSegment<double> lsa{geo::Point<double>(1, 1),
                               geo::Point<double>(2, 2)};
  geo::LineSegment<double> lsb{geo::Point<double>(1, 2),
                               geo::Point<double>(2, 1)};
  geo::LineSegment<double> lsc{geo::Point<double>(2.1, 2),
                               geo::Point<double>(3, 3)};

  TEST(geo::crossProd(lsa.first, lsb), ==, approx(-1));
  TEST(geo::crossProd(lsa.second, lsb), ==, approx(1));

  TEST(geo::intersects(lsa, lsb));

  TEST(!geo::intersects(lsa, lsa));
  TEST(!geo::intersects(lsb, lsb));
  TEST(!geo::intersects(lsa, lsc));

  TEST(!geo::intersects(geo::Point<double>(871569.2, 6104550.4),
                          geo::Point<double>(871581.2, 6104536),
                          geo::Point<double>(871580.3, 6104541.3),
                          geo::Point<double>(871625.7, 6104510.1)));

  TEST(!geo::intersects(geo::Point<double>(0, 0), geo::Point<double>(1, 1),
                          geo::Point<double>(0.5, 0.5),
                          geo::Point<double>(1.5, 1.5)));

  geo::Line<double> l{geo::Point<double>(1, 1), geo::Point<double>(2, 2),
                      geo::Point<double>(2, 4)};
  TEST(!geo::contains(geo::Point<double>(1, 2), l));
  TEST(geo::contains(geo::Point<double>(2, 2), l));
  TEST(geo::contains(geo::Point<double>(2, 3), l));

  geo::Box<double> bbox(geo::Point<double>(1, 1), geo::Point<double>(3, 3));
  TEST(geo::intersects(l, bbox));
  geo::Line<double> ll{geo::Point<double>(0, 0), geo::Point<double>(4, 4)};
  TEST(geo::intersects(ll, bbox));
  geo::Line<double> lll{geo::Point<double>(0, 0), geo::Point<double>(0, 4)};
  TEST(!geo::intersects(lll, bbox));
  geo::Line<double> llll{geo::Point<double>(1.2, 0), geo::Point<double>(1, 2)};
  TEST(geo::intersects(llll, bbox));

  Line<double> l5new;
  l5new.push_back(Point<double>(-10, -5));
  l5new.push_back(Point<double>(-8, -4));
  TEST(geo::getBoundingBox(l5new).getUpperRight().getX(), ==, approx(-8));
  TEST(geo::getBoundingBox(l5new).getUpperRight().getY(), ==, approx(-4));

  Line<double> l5;
  l5.push_back(Point<double>(0, 0));
  l5.push_back(Point<double>(1.5, 2));
  Box<double> req(Point<double>(.5, 1), Point<double>(1, 1.5));

  TEST(geo::getBoundingBox(l5[0]).getLowerLeft().getX(), ==, approx(0));
  TEST(geo::getBoundingBox(l5[0]).getLowerLeft().getY(), ==, approx(0));

  TEST(geo::getBoundingBox(l5).getLowerLeft().getX(), ==, approx(0));
  TEST(geo::getBoundingBox(l5).getLowerLeft().getY(), ==, approx(0));
  TEST(geo::getBoundingBox(l5).getUpperRight().getX(), ==, approx(1.5));
  TEST(geo::getBoundingBox(l5).getUpperRight().getY(), ==, approx(2));
  TEST(geo::intersects(geo::getBoundingBox(l5),
                         geo::getBoundingBox(Line<double>{
                             Point<double>(.5, 1), Point<double>(1, 1)})));
  TEST(geo::intersects(
      l5, Line<double>{Point<double>(.5, 1), Point<double>(1, 1)}));
  TEST(geo::intersects(l5, req));

  Box<double> boxa(Point<double>(1, 1), Point<double>(2, 2));
  TEST(geo::intersects(
      boxa, Box<double>(Point<double>(1.5, 1.5), Point<double>(1.7, 1.7))));
  TEST(geo::intersects(
      boxa, Box<double>(Point<double>(0, 0), Point<double>(3, 3))));
  TEST(geo::intersects(
      boxa, Box<double>(Point<double>(1.5, 1.5), Point<double>(3, 3))));
  TEST(geo::intersects(
      boxa, Box<double>(Point<double>(0, 0), Point<double>(1.5, 1.5))));

  TEST(geo::intersects(
      Box<double>(Point<double>(1.5, 1.5), Point<double>(1.7, 1.7)), boxa));
  TEST(geo::intersects(Box<double>(Point<double>(0, 0), Point<double>(3, 3)),
                         boxa));
  TEST(geo::intersects(
      Box<double>(Point<double>(1.5, 1.5), Point<double>(3, 3)), boxa));
  TEST(geo::intersects(
      Box<double>(Point<double>(0, 0), Point<double>(1.5, 1.5)), boxa));

  Polygon<double> poly({Point<double>(1, 1), Point<double>(3, 2),
                        Point<double>(4, 3), Point<double>(6, 3),
                        Point<double>(5, 1)});
  TEST(geo::getWKT(poly), ==, "POLYGON ((1 1, 3 2, 4 3, 6 3, 5 1, 1 1))");
  TEST(geo::contains(Point<double>(4, 2), poly));
  TEST(!geo::contains(Point<double>(3, 3), poly));
  TEST(geo::contains(Point<double>(1, 1), poly));
  TEST(geo::contains(Point<double>(3, 2), poly));
  TEST(geo::contains(Point<double>(4, 3), poly));
  TEST(geo::contains(Point<double>(6, 3), poly));
  TEST(geo::contains(Point<double>(5, 1), poly));

  TEST(geo::contains(Line<double>{Point<double>(6, 3), Point<double>(5, 1)},
                       poly));
  TEST(!geo::contains(Line<double>{Point<double>(6, 3), Point<double>(50, 1)},
                        poly));
  TEST(geo::contains(Line<double>{Point<double>(4, 2), Point<double>(4.5, 2)},
                       poly));
  TEST(geo::contains(Line<double>{Point<double>(4, 2), Point<double>(5, 1)},
                       poly));

  Box<double> polybox(Point<double>(1, 1), Point<double>(6, 4));
  TEST(geo::centroid(polybox).getX(), ==, approx(3.5));
  TEST(geo::centroid(polybox).getY(), ==, approx(2.5));
  TEST(geo::contains(poly, polybox));
  TEST(!geo::contains(polybox, poly));
  Box<double> polybox2(Point<double>(4, 1), Point<double>(5, 2));
  TEST(geo::contains(polybox2, poly));
  TEST(geo::contains(poly, getBoundingBox(poly)));

  Point<double> rotP(2, 2);
  TEST(geo::dist(geo::rotate(rotP, 180, Point<double>(1, 1)),
                   Point<double>(0, 0)) == approx(0));
  TEST(geo::dist(geo::rotate(rotP, 360, Point<double>(1, 1)), rotP) ==
         approx(0));

  Line<double> rotLine({{1, 1}, {3, 3}});
  TEST(geo::rotate(rotLine, 90, Point<double>(2, 2))[0].getX(), ==, approx(1));
  TEST(geo::rotate(rotLine, 90, Point<double>(2, 2))[0].getY(), ==, approx(3));
  TEST(geo::rotate(rotLine, 90, Point<double>(2, 2))[1].getX(), ==, approx(3));
  TEST(geo::rotate(rotLine, 90, Point<double>(2, 2))[1].getY(), ==, approx(1));

  MultiLine<double> multiRotLine({{{1, 1}, {3, 3}}, {{1, 3}, {3, 1}}});
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][0].getX() ==
         approx(1));
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][0].getY() ==
         approx(3));
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][1].getX() ==
         approx(3));
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][1].getY() ==
         approx(1));
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][0].getX() ==
         approx(3));
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][0].getY() ==
         approx(3));
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][1].getX() ==
         approx(1));
  TEST(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][1].getY() ==
         approx(1));

  TEST(geo::getWKT(multiRotLine) ==
         "MULTILINESTRING ((1 1, 3 3), (1 3, 3 1))");

  TEST(geo::contains(
      multiRotLine[0],
      geo::move(geo::move(multiRotLine, 1.0, 2.0), -1.0, -2.0)[0]));
  TEST(geo::contains(multiRotLine, geo::getBoundingBox(Line<double>{
                                         {1, 1}, {3, 3}, {1, 3}, {3, 1}})));

  TEST(geo::contains(
      getBoundingBox(multiRotLine),
      geo::getBoundingBox(Line<double>{{1, 1}, {3, 3}, {1, 3}, {3, 1}})));
  TEST(geo::contains(
      geo::getBoundingBox(Line<double>{{1, 1}, {3, 3}, {1, 3}, {3, 1}}),
      getBoundingBox(multiRotLine)));

  TEST(geo::dist(geo::centroid(rotP), rotP), ==, approx(0));
  TEST(geo::dist(geo::centroid(rotLine), rotP), ==, approx(0));
  TEST(geo::dist(geo::centroid(polybox), Point<double>(3.5, 2.5)) ==
         approx(0));
  TEST(geo::dist(geo::centroid(Polygon<double>({{0, 0}, {3, 4}, {4, 3}})),
                   Point<double>(7.0 / 3.0, 7.0 / 3.0)) == approx(0));

  auto polyy = Polygon<double>({{0, 0}, {3, 4}, {4, 3}});
  MultiPolygon<double> mpoly{polyy, polyy};

  TEST(geo::getWKT(polyy), ==, "POLYGON ((0 0, 3 4, 4 3, 0 0))");
  TEST(geo::getWKT(mpoly) ==
         "MULTIPOLYGON (((0 0, 3 4, 4 3, 0 0)), ((0 0, 3 4, 4 3, 0 0)))");

  auto hull = geo::convexHull(Line<double>{
      {0.1, 3}, {1, 1}, {2, 2}, {4, 4}, {0, 0}, {1, 2}, {3, 1}, {3, 3}});
  TEST(hull.getOuter().size(), ==, size_t(4));
  TEST(hull.getOuter()[0].getX(), ==, approx(0));
  TEST(hull.getOuter()[0].getY(), ==, approx(0));
  TEST(hull.getOuter()[1].getX(), ==, approx(0.1));
  TEST(hull.getOuter()[1].getY(), ==, approx(3));
  TEST(hull.getOuter()[2].getX(), ==, approx(4));
  TEST(hull.getOuter()[2].getY(), ==, approx(4));
  TEST(hull.getOuter()[3].getX(), ==, approx(3));
  TEST(hull.getOuter()[3].getY(), ==, approx(1));
  TEST(geo::contains(geo::convexHull(geo::getBoundingBox(poly)),
                       geo::getBoundingBox(poly)));
  TEST(geo::contains(geo::getBoundingBox(poly),
                       geo::convexHull(geo::getBoundingBox(poly))));

  auto hull2 = geo::convexHull(Line<double>{{0.1, 3},
                                            {1, 1},
                                            {2, 2},
                                            {4, 4},
                                            {0, 0},
                                            {1, 2},
                                            {3, 1},
                                            {3, 3},
                                            {-0.1, 1}});
  TEST(hull2.getOuter().size(), ==, size_t(5));
  TEST(hull2.getOuter()[0].getX(), ==, approx(-.1));
  TEST(hull2.getOuter()[0].getY(), ==, approx(1));
  TEST(hull2.getOuter()[1].getX(), ==, approx(0.1));
  TEST(hull2.getOuter()[1].getY(), ==, approx(3));
  TEST(hull2.getOuter()[2].getX(), ==, approx(4));
  TEST(hull2.getOuter()[2].getY(), ==, approx(4));
  TEST(hull2.getOuter()[3].getX(), ==, approx(3));
  TEST(hull2.getOuter()[3].getY(), ==, approx(1));
  TEST(hull2.getOuter()[4].getX(), ==, approx(0));
  TEST(hull2.getOuter()[4].getY(), ==, approx(0));

  auto hull3 =
      geo::convexHull(Line<double>{{0.1, 3}, {4, 4}, {0, 0}, {1, 2}, {3, 1}});
  TEST(hull3.getOuter().size(), ==, size_t(4));
  TEST(hull3.getOuter()[0].getX(), ==, approx(0));
  TEST(hull3.getOuter()[0].getY(), ==, approx(0));
  TEST(hull3.getOuter()[3].getX(), ==, approx(3));
  TEST(hull3.getOuter()[3].getY(), ==, approx(1));
  TEST(hull3.getOuter()[2].getX(), ==, approx(4));
  TEST(hull3.getOuter()[2].getY(), ==, approx(4));
  TEST(hull3.getOuter()[1].getX(), ==, approx(0.1));
  TEST(hull3.getOuter()[1].getY(), ==, approx(3));

  hull3 = geo::convexHull(
      Line<double>{{0.1, 3}, {4, 4}, {2, 1}, {3, 2}, {0, 0}, {1, 2}, {3, 1}});
  TEST(hull3.getOuter().size(), ==, size_t(4));
  TEST(hull3.getOuter()[0].getX(), ==, approx(0));
  TEST(hull3.getOuter()[0].getY(), ==, approx(0));
  TEST(hull3.getOuter()[3].getX(), ==, approx(3));
  TEST(hull3.getOuter()[3].getY(), ==, approx(1));
  TEST(hull3.getOuter()[2].getX(), ==, approx(4));
  TEST(hull3.getOuter()[2].getY(), ==, approx(4));
  TEST(hull3.getOuter()[1].getX(), ==, approx(0.1));
  TEST(hull3.getOuter()[1].getY(), ==, approx(3));

  hull3 = geo::convexHull(Line<double>{
      {4, 4}, {1, 2}, {2, 1}, {3, 2}, {0.1, 3}, {0, 0}, {1, 2}, {3, 1}});
  TEST(hull3.getOuter().size(), ==, size_t(4));
  TEST(hull3.getOuter()[0].getX(), ==, approx(0));
  TEST(hull3.getOuter()[0].getY(), ==, approx(0));
  TEST(hull3.getOuter()[3].getX(), ==, approx(3));
  TEST(hull3.getOuter()[3].getY(), ==, approx(1));
  TEST(hull3.getOuter()[2].getX(), ==, approx(4));
  TEST(hull3.getOuter()[2].getY(), ==, approx(4));
  TEST(hull3.getOuter()[1].getX(), ==, approx(0.1));
  TEST(hull3.getOuter()[1].getY(), ==, approx(3));

  hull3 = geo::convexHull(Line<double>{{4, 4}, {1, 2}, {3, 1}});
  TEST(hull3.getOuter().size(), ==, size_t(3));
  TEST(hull3.getOuter()[0].getX(), ==, approx(1));
  TEST(hull3.getOuter()[0].getY(), ==, approx(2));
  TEST(hull3.getOuter()[2].getX(), ==, approx(3));
  TEST(hull3.getOuter()[2].getY(), ==, approx(1));
  TEST(hull3.getOuter()[1].getX(), ==, approx(4));
  TEST(hull3.getOuter()[1].getY(), ==, approx(4));

  hull3 = geo::convexHull(Line<double>{{4, 4}, {1, 2}, {3, 10}});
  TEST(hull3.getOuter().size(), ==, size_t(3));
  TEST(hull3.getOuter()[0].getX(), ==, approx(1));
  TEST(hull3.getOuter()[0].getY(), ==, approx(2));
  TEST(hull3.getOuter()[2].getX(), ==, approx(4));
  TEST(hull3.getOuter()[2].getY(), ==, approx(4));
  TEST(hull3.getOuter()[1].getX(), ==, approx(3));
  TEST(hull3.getOuter()[1].getY(), ==, approx(10));

  Line<double> test{{0.3215348546593775, 0.03629583077160248},
                    {0.02402358131857918, -0.2356728797179394},
                    {0.04590851212470659, -0.4156409924995536},
                    {0.3218384001607433, 0.1379850698988746},
                    {0.11506479756447, -0.1059521474930943},
                    {0.2622539999543261, -0.29702873322836},
                    {-0.161920957418085, -0.4055339716426413},
                    {0.1905378631228002, 0.3698601009043493},
                    {0.2387090918968516, -0.01629827079949742},
                    {0.07495888748668034, -0.1659825110491202},
                    {0.3319341836794598, -0.1821814101954749},
                    {0.07703635755650362, -0.2499430638271785},
                    {0.2069242999022122, -0.2232970760420869},
                    {0.04604079532068295, -0.1923573186549892},
                    {0.05054295812784038, 0.4754929463150845},
                    {-0.3900589168910486, 0.2797829520700341},
                    {0.3120693385713448, -0.0506329867529059},
                    {0.01138812723698857, 0.4002504701728471},
                    {0.009645149586391732, 0.1060251100976254},
                    {-0.03597933197019559, 0.2953639456959105},
                    {0.1818290866742182, 0.001454397571696298},
                    {0.444056063372694, 0.2502497166863175},
                    {-0.05301752458607545, -0.06553921621808712},
                    {0.4823896228171788, -0.4776170002088109},
                    {-0.3089226845734964, -0.06356112199235814},
                    {-0.271780741188471, 0.1810810595574612},
                    {0.4293626522918815, 0.2980897964891882},
                    {-0.004796652127799228, 0.382663812844701},
                    {0.430695573269106, -0.2995073500084759},
                    {0.1799668387323309, -0.2973467472915973},
                    {0.4932166845474547, 0.4928094162538735},
                    {-0.3521487911717489, 0.4352656197131292},
                    {-0.4907368011686362, 0.1865826865533206},
                    {-0.1047924716070224, -0.247073392148198},
                    {0.4374961861758457, -0.001606279519951237},
                    {0.003256207800708899, -0.2729194320486108},
                    {0.04310378203457577, 0.4452604050238248},
                    {0.4916198379282093, -0.345391701297268},
                    {0.001675087028811806, 0.1531837672490476},
                    {-0.4404289572876217, -0.2894855991839297}

  };
  hull3 = geo::convexHull(test);
  TEST(geo::contains(test, hull3));
  TEST(hull3.getOuter().size(), ==, size_t(8));
  TEST(geo::contains(
      Polygon<double>({{-0.161920957418085, -0.4055339716426413},
                       {0.05054295812784038, 0.4754929463150845},
                       {0.4823896228171788, -0.4776170002088109},
                       {0.4932166845474547, 0.4928094162538735},
                       {-0.3521487911717489, 0.4352656197131292},
                       {-0.4907368011686362, 0.1865826865533206},
                       {0.4916198379282093, -0.345391701297268},
                       {-0.4404289572876217, -0.2894855991839297}}),
      hull3));
  TEST(geo::contains(
      hull3, Polygon<double>({{-0.161920957418085, -0.4055339716426413},
                              {0.05054295812784038, 0.4754929463150845},
                              {0.4823896228171788, -0.4776170002088109},
                              {0.4932166845474547, 0.4928094162538735},
                              {-0.3521487911717489, 0.4352656197131292},
                              {-0.4907368011686362, 0.1865826865533206},
                              {0.4916198379282093, -0.345391701297268},
                              {-0.4404289572876217, -0.2894855991839297}})));

  hull3 = geo::convexHull(Line<double>{{3, 6},
                                       {8, 10},
                                       {3, 5},
                                       {20, -10},
                                       {-4, 5},
                                       {10, 2},
                                       {5, 1},
                                       {45, 1},
                                       {30, -9},
                                       {3, 14},
                                       {25, -5.5}});
  TEST(hull3.getOuter().size(), ==, size_t(5));
  TEST(hull3.getOuter()[0].getX(), ==, approx(-4));
  TEST(hull3.getOuter()[0].getY(), ==, approx(5));
  TEST(hull3.getOuter()[4].getX(), ==, approx(20));
  TEST(hull3.getOuter()[4].getY(), ==, approx(-10));
  TEST(hull3.getOuter()[3].getX(), ==, approx(30));
  TEST(hull3.getOuter()[3].getY(), ==, approx(-9));
  TEST(hull3.getOuter()[2].getX(), ==, approx(45));
  TEST(hull3.getOuter()[2].getY(), ==, approx(1));
  TEST(hull3.getOuter()[1].getX(), ==, approx(3));
  TEST(hull3.getOuter()[1].getY(), ==, approx(14));

  hull3 = geo::convexHull(Line<double>{
      {7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}});
  TEST(hull3.getOuter().size(), ==, size_t(8));
  TEST(geo::contains(geo::Polygon<double>({{-9, 0},
                                             {-7, -7},
                                             {0, -9},
                                             {7, -7},
                                             {9, 0},
                                             {7, 7},
                                             {0, 9},
                                             {-7, 7}}),
                       hull3));
  TEST(geo::contains(hull3, geo::Polygon<double>({{-9, 0},
                                                    {-7, -7},
                                                    {0, -9},
                                                    {7, -7},
                                                    {9, 0},
                                                    {7, 7},
                                                    {0, 9},
                                                    {-7, 7}})));

  hull3 = geo::convexHull(Line<double>{{7, 7},
                                       {7, -7},
                                       {-7, -7},
                                       {-7, 7},
                                       {9, 0},
                                       {-9, 0},
                                       {0, 9},
                                       {0, -9},
                                       {0, 0},
                                       {1, 2},
                                       {-2, 1},
                                       {-1, -1},
                                       {3, 4},
                                       {4, 3},
                                       {-5, 4},
                                       {6, 5}});
  TEST(hull3.getOuter().size(), ==, size_t(8));
  TEST(geo::contains(geo::Polygon<double>({{-9, 0},
                                             {-7, -7},
                                             {0, -9},
                                             {7, -7},
                                             {9, 0},
                                             {7, 7},
                                             {0, 9},
                                             {-7, 7}}),
                       hull3));
  TEST(geo::contains(hull3, geo::Polygon<double>({{-9, 0},
                                                    {-7, -7},
                                                    {0, -9},
                                                    {7, -7},
                                                    {9, 0},
                                                    {7, 7},
                                                    {0, 9},
                                                    {-7, 7}})));

  hull3 = geo::convexHull(Line<double>{
      {0, 0},   {1, 2},  {-2, 1}, {-1, -1}, {3, 4},   {4, 3},   {-5, 4},
      {6, 5},   {7, 7},  {7, -7}, {-7, -7}, {-7, 7},  {9, 0},   {-9, 0},
      {0, 9},   {0, -9}, {-8, 0}, {8, 0},   {-7, 0},  {7, 0},   {-6, 0},
      {6, 0},   {-5, 0}, {5, 0},  {-4, 0},  {4, 0},   {-3, 0},  {3, 0},
      {-2, 0},  {2, 0},  {-1, 0}, {1, 0},   {0, -8},  {0, 8},   {0, -7},
      {0, 7},   {0, -6}, {0, 6},  {0, -5},  {0, 5},   {0, -4},  {0, 4},
      {0, -3},  {0, 3},  {0, -2}, {0, 2},   {0, -1},  {0, 1},   {1, 1},
      {2, 2},   {3, 3},  {4, 4},  {5, 5},   {6, 6},   {1, -1},  {2, -2},
      {3, -3},  {4, -4}, {5, -5}, {6, -6},  {-1, 1},  {-2, 2},  {-3, 3},
      {-4, 4},  {-5, 5}, {-6, 6}, {-1, -1}, {-2, -2}, {-3, -3}, {-4, -4},
      {-5, -5}, {-6, -6}});
  TEST(hull3.getOuter().size(), ==, size_t(8));
  TEST(geo::contains(geo::Polygon<double>({{-9, 0},
                                             {-7, -7},
                                             {0, -9},
                                             {7, -7},
                                             {9, 0},
                                             {7, 7},
                                             {0, 9},
                                             {-7, 7}}),
                       hull3));
  TEST(geo::contains(hull3, geo::Polygon<double>({{-9, 0},
                                                    {-7, -7},
                                                    {0, -9},
                                                    {7, -7},
                                                    {9, 0},
                                                    {7, 7},
                                                    {0, 9},
                                                    {-7, 7}})));

  TEST(geo::area(geo::Point<double>(1, 2)), ==, approx(0));
  TEST(geo::area(geo::Line<double>{{1, 2}, {2, 5}}), ==, approx(0));
  TEST(geo::area(geo::Box<double>({0, 0}, {1, 1})), ==, approx(1));
  TEST(geo::area(geo::Box<double>({1, 1}, {1, 1})), ==, approx(0));
  TEST(geo::area(geo::Box<double>({0, 0}, {2, 2})), ==, approx(4));
  TEST(geo::area(geo::Polygon<double>({{0, 0}, {1, 0}, {1, 1}, {0, 1}})) ==
         approx(1));
  TEST(geo::area(geo::Polygon<double>({{0, 0}, {1, 0}, {1, 1}})) ==
         approx(0.5));

  auto obox =
      geo::getOrientedEnvelope(geo::Line<double>{{0, 0}, {1, 1}, {1.5, 0.5}});
  TEST(geo::contains(
      geo::convexHull(obox),
      geo::Polygon<double>({{0.0, 0.0}, {1.0, 1.0}, {1.5, 0.5}, {0.5, -0.5}})));
  TEST(geo::contains(
      geo::Polygon<double>({{0.0, 0.0}, {1.0, 1.0}, {1.5, 0.5}, {0.5, -0.5}}),
      geo::convexHull(obox)));

  TEST(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}},
                   geo::LineSegment<double>{{2, 2}, {2, 0}}) == approx(0));
  TEST(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}},
                   geo::LineSegment<double>{{2, 4}, {2, 2}}) == approx(1));
  TEST(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}},
                   geo::LineSegment<double>{{1, 1}, {3, 1}}) == approx(0));
  TEST(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}},
                   geo::LineSegment<double>{{1, 2}, {3, 2}}) == approx(1));
  TEST(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}},
                   geo::LineSegment<double>{{1, 2}, {3, 5}}) == approx(1));

  TEST(geo::dist(geo::Line<double>{{1, 1}, {3, 1}},
                   geo::Point<double>{2, 1}) == approx(0));
  TEST(geo::dist(geo::Line<double>{{1, 1}, {3, 1}},
                   geo::Point<double>{2, 2}) == approx(1));
  TEST(geo::dist(geo::Line<double>{{1, 1}, {3, 1}},
                   geo::Point<double>{3, 1}) == approx(0));
  TEST(geo::dist(geo::Line<double>{{1, 1}, {3, 1}},
                   geo::Point<double>{1, 1}) == approx(0));

  TEST(geo::dist(Line<double>{{7, 7},
                                {7, -7},
                                {-7, -7},
                                {-7, 7},
                                {9, 0},
                                {-9, 0},
                                {0, 9},
                                {0, -9}},
                   Line<double>{{7, 7},
                                {7, -7},
                                {-7, -7},
                                {-7, 7},
                                {9, 0},
                                {-9, 0},
                                {0, 9},
                                {0, -9}}) == approx(0));
  TEST(geo::dist(Line<double>{{7, 7},
                                {7, -7},
                                {-7, -7},
                                {-7, 7},
                                {9, 0},
                                {-9, 0},
                                {0, 9},
                                {0, -9}},
                   LineSegment<double>{{6, 7}, {8, -7}}) == approx(0));
  TEST(geo::dist(Line<double>{{7, 7},
                                {7, -7},
                                {-7, -7},
                                {-7, 7},
                                {9, 0},
                                {-9, 0},
                                {0, 9},
                                {0, -9}},
                   Point<double>{7, 4}) == approx(0));
  TEST(geo::dist(Line<double>{{0, 0}, {1, 1}, {2, 0}},
                   Line<double>{{1.5, 0.5}, {1.5, 100}}) == approx(0));
  TEST(geo::dist(Line<double>{{0, 0}, {1, 1}, {2, 0}},
                   Line<double>{{2, 0.5}, {2, 100}}) == approx(0.353553));

  TEST(geo::contains(util::geo::Point<double>{1.5, 0.5},
                       util::geo::LineSegment<double>{{1, 1}, {1.5, 0.5}}));
  TEST(geo::contains(util::geo::Point<double>{1.5, 0.5},
                       util::geo::LineSegment<double>{{1, 1}, {1.5, 0.5}}));

  auto polyTest =
      geo::Polygon<double>({{1, 1}, {3, 1}, {2, 2}, {3, 3}, {1, 3}});
  TEST(!geo::contains(util::geo::LineSegment<double>({2.5, 1.3}, {2.5, 2.6}),
                        polyTest));

  TEST(!geo::contains(util::geo::LineSegment<double>{{2.5, 1.3}, {2.5, 2.6}},
                        polyTest));
  TEST(geo::contains(util::geo::LineSegment<double>{{2.5, 2.6}, {1.5, 2}},
                       polyTest));
  TEST(!geo::contains(
      util::geo::Line<double>{{2.5, 1.3}, {2.5, 2.6}, {1.5, 2}}, polyTest));
  TEST(geo::contains(
      util::geo::Line<double>{{2.5, 1.3}, {1.5, 2}, {2.5, 2.6}}, polyTest));

  TEST(!geo::contains(util::geo::Box<double>{{1, 1}, {2.5, 2.6}}, polyTest));

  TEST(
      geo::intersects(Box<double>(Point<double>(0, 0), Point<double>(10, 10)),
                      Box<double>(Point<double>(2, 2), Point<double>(8, 8))));
  TEST(
      geo::intersects(Box<double>(Point<double>(0, 0), Point<double>(10, 10)),
                      Box<double>(Point<double>(-2, -2), Point<double>(8, 8))));
  TEST(geo::intersects(
      Box<double>(Point<double>(0, 0), Point<double>(10, 10)),
      Box<double>(Point<double>(-2, -2), Point<double>(12, 12))));
  TEST(
      geo::intersects(Box<double>(Point<double>(0, 0), Point<double>(10, 10)),
                      Box<double>(Point<double>(5, 5), Point<double>(12, 12))));

  TEST(!geo::intersects(
      Box<double>(Point<double>(0, 0), Point<double>(10, 10)),
      Box<double>(Point<double>(15, 15), Point<double>(12, 12))));

  double rad = 10.0;
  int n = 20;
  util::geo::MultiPoint<double> mp;

  for (int i = 0; i < n; i++) {
    double x = rad * cos((2.0 * M_PI / static_cast<double>(n)) *
                         static_cast<double>(i));
    double y = rad * sin((2.0 * M_PI / static_cast<double>(n)) *
                         static_cast<double>(i));

    mp.push_back(util::geo::DPoint(x, y));
  }

  auto h = util::geo::convexHull(mp);

  TEST(geo::contains(mp, h));

  TEST(geo::dist(geo::interpolate(DPoint(1.0, 1.0), DPoint(1.0, 3.0), 0.5), DPoint(1.0, 2.0)), ==, approx(0));
  TEST(geo::dist(geo::interpolate(DPoint(1.0, 1.0), DPoint(2.0, 2.0), 0), DPoint(1.0, 1.0)), ==, approx(0));
  TEST(geo::dist(geo::interpolate(DPoint(1.0, 1.0), DPoint(2.0, 2.0), 0.5), DPoint(1.5, 1.5)), ==, approx(0));
  TEST(geo::dist(geo::interpolate(DPoint(1.0, 1.0), DPoint(2.0, 2.0), 1), DPoint(2, 2)), ==, approx(0));
  TEST(geo::dist(geo::pointAtDist(DLine{{0, 0}, {0, 1}, {0, 2}}, 1), DPoint{0, 1}), ==, approx(0));
  TEST(geo::dist(geo::pointAtDist(DLine{{0, 0}, {0, 1}, {0, 2}}, 2), DPoint{0, 2}), ==, approx(0));
  TEST(geo::dist(geo::pointAtDist(DLine{{0, 0}, {0, 1}, {0, 2}}, 0), DPoint{0, 0}), ==, approx(0));

  TEST(geo::getWKT(geo::orthoLineAtDist(DLine{{0, 0}, {0, 1}, {0, 2}}, 1, 1)), ==, "LINESTRING (-0.5 1, 0.5 1)");

  TEST(geo::getWKT(geo::segment(DLine{{0, 0}, {0, 1}, {0, 2}}, 0, 0.5)), ==, "LINESTRING (0 0, 0 1)");
  TEST(geo::getWKT(geo::segment(DLine{{0, 0}, {0, 1}, {0, 2}}, 0.5, 1)), ==, "LINESTRING (0 1, 0 2)");
}

  // inversion count
  std::vector<int> test = {2, 1};
  TEST(inversions(test), ==, 1);

  test = {};
  TEST(inversions(test), ==, 0);

  test = {2};
  TEST(inversions(test), ==, 0);

  test = {2, 1};
  TEST(inversions(test), ==, 1);

  test = {1, 2};
  TEST(inversions(test), ==, 0);

  test = {2, 1, 3};
  TEST(inversions(test), ==, 1);

  test = {2, 3, 1};
  TEST(inversions(test), ==, 2);

  test = {3, 2, 1};
  TEST(inversions(test), ==, 3);

  test = {1, 2, 3};
  TEST(inversions(test), ==, 0);

  test = {1, 3, 2, 6, 5, 4, 8, 7, 9};
  TEST(inversions(test), ==, 5);

  test = {1, 2, 3, 4, 5, 6};
  TEST(inversions(test), ==, 0);

  test = {9, 8, 7, 6, 5, 4, 3, 2, 1};
  TEST(inversions(test), ==, 8 + 7 + 6 + 5 + 4 + 3 + 2 + 1);

  // nice float formatting
	TEST(formatFloat(15.564, 3), ==, "15.564");
	TEST(formatFloat(15.564, 0), ==, "16");
	TEST(formatFloat(15.0000, 10), ==, "15");
	TEST(formatFloat(15.0100, 10), ==, "15.01");
	TEST(formatFloat(0.0000, 10), ==, "0");
	TEST(formatFloat(-1.0000, 10), ==, "-1");
	TEST(formatFloat(-15.000001, 10), ==, "-15.000001");
}
