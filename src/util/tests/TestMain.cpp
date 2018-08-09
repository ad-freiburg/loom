// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include "lest.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/json/Writer.h"
#include "util/graph/DirGraph.h"
#include "util/graph/UndirGraph.h"
#include "util/graph/Dijkstra.h"
#include "util/graph/EDijkstra.h"
#include "util/geo/Grid.h"
#include "util/Misc.h"

using lest::approx;
using namespace util;
using namespace util::geo;
using namespace util::graph;

// define LEST cases
const lest::test specification[] = {

// ___________________________________________________________________________
{
CASE("atof") {
  EXPECT(util::atof("45.534215") == approx(45.534215));
  EXPECT(util::atof("5.534") == approx(5.534));
  EXPECT(util::atof("534") == approx(534));
  EXPECT(util::atof("-534") == approx(-534));
  EXPECT(util::atof("-45.534215") == approx(-45.534215));
  EXPECT(util::atof("-45.534215", 2) == approx(-45.53));


  // TODO: more test cases
}},

// ___________________________________________________________________________
{
CASE("json") {
  std::stringstream ss;
  util::json::Writer wr(&ss, 2, false);

  util::json::Val a("bla");
  util::json::Val b(1);
  util::json::Val c(1.0);
  util::json::Val d("a");
  util::json::Val e({"a", "b", "c"});

  util::json::Val f({1, json::Array{2, 3, 4}, 3});

  ss = std::stringstream();
  wr = util::json::Writer(&ss, 2, false);
  util::json::Val i({1, json::Array{2, json::Null(), 4}, true});
  wr.val(i);
  wr.closeAll();
  EXPECT(ss.str() == "[1,[2,null,4],true]");

  ss = std::stringstream();
  wr = util::json::Writer(&ss, 2, false);
  i = util::json::Val({1, json::Array{2, json::Null(), 4}, false});
  wr.val(i);
  wr.closeAll();
  EXPECT(ss.str() == "[1,[2,null,4],false]");

  ss = std::stringstream();
  wr = util::json::Writer(&ss, 2, false);
  i = util::json::Val({1, json::Array{2, json::Null(), 4}, false});
  wr.val(i);
  wr.closeAll();
  EXPECT(ss.str() == "[1,[2,null,4],false]");

  ss = std::stringstream();
  wr = util::json::Writer(&ss, 2, false);
  i = util::json::Val({1, json::Array{2.13, "", 4}, 0});
  wr.val(i);
  wr.closeAll();
  EXPECT(ss.str() == "[1,[2.13,\"\",4],0]");

  ss = std::stringstream();
  wr = util::json::Writer(&ss, 2, false);
  i = util::json::Val({1, json::Array{2.13, json::Dict{{"a", 1}, {"B", 2.123}}, 4}, 0});
  wr.val(i);
  wr.closeAll();
  EXPECT((ss.str() == "[1,[2.13,{\"a\":1,\"B\":2.12},4],0]" || ss.str() == "[1,[2.13,{\"B\":2.12,\"a\":1},4],0]"));
}},

// ___________________________________________________________________________
{
CASE("dirgraph") {
  DirGraph<int, int> g;

  DirNode<int, int>* a = new DirNode<int, int>(0);
  DirNode<int, int>* b = new DirNode<int, int>(0);
  g.addNd(a);
  EXPECT(g.getNds()->size() == (size_t)1);
  g.addNd(b);
  EXPECT(g.getNds()->size() == (size_t)2);

  g.addEdg(a, b);
  EXPECT(a->getDeg() == (size_t)1);
  EXPECT(b->getDeg() == (size_t)0);

  auto c = g.addNd();

  g.addEdg(a, c);
  g.addEdg(c, b);
  EXPECT(a->getDeg() == (size_t)2);
  EXPECT(b->getDeg() == (size_t)0);
  EXPECT(c->getDeg() == (size_t)1);

  g.delEdg(a, c);

  EXPECT(a->getDeg() == (size_t)1);
  EXPECT(b->getDeg() == (size_t)0);
  EXPECT(c->getDeg() == (size_t)1);

  g.addEdg(a, a);
  EXPECT(a->getDeg() == (size_t)2);

  g.delEdg(a, a);
  EXPECT(a->getDeg() == (size_t)1);

  g.delEdg(a, a);
  EXPECT(a->getDeg() == (size_t)1);

  // TODO: more test cases
}},

// ___________________________________________________________________________
{
CASE("unddirgraph") {
  UndirGraph<int, int> g;

  UndirNode<int, int>* a = new UndirNode<int, int>(0);
  UndirNode<int, int>* b = new UndirNode<int, int>(0);
  g.addNd(a);
  EXPECT(g.getNds()->size() == (size_t)1);
  g.addNd(b);
  EXPECT(g.getNds()->size() == (size_t)2);

  g.addEdg(a, b);
  EXPECT(a->getDeg() == (size_t)1);
  EXPECT(b->getDeg() == (size_t)1);

  auto c = g.addNd();

  g.addEdg(a, c);
  g.addEdg(c, b);
  EXPECT(a->getDeg() == (size_t)2);
  EXPECT(b->getDeg() == (size_t)2);
  EXPECT(c->getDeg() == (size_t)2);

  g.delEdg(a, c);

  EXPECT(a->getDeg() == (size_t)1);
  EXPECT(b->getDeg() == (size_t)2);
  EXPECT(c->getDeg() == (size_t)1);

  g.delNd(b);

  EXPECT(a->getDeg() == (size_t)0);
  EXPECT(c->getDeg() == (size_t)0);

  g.addEdg(a, a);
  EXPECT(a->getDeg() == (size_t)1);

  g.delEdg(a, a);
  EXPECT(a->getDeg() == (size_t)0);

  // TODO: more test cases

}},

// ___________________________________________________________________________
{
CASE("grid") {
  Grid<int, Line, double> g(.5, .5, Box<double>(Point<double>(0, 0), Point<double>(3, 3)));

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
  EXPECT(ret.size() == (size_t)1);

  ret.clear();
  g.getNeighbors(1, 0, &ret);
  EXPECT(ret.size() == (size_t)1);

  ret.clear();
  g.getNeighbors(1, 0.55, &ret);
  EXPECT(ret.size() == (size_t)2);


  // TODO: more test cases
}},

// ___________________________________________________________________________
{
CASE("densify") {
  Line<double> a;
  a.push_back(Point<double>(1, 1));
  a.push_back(Point<double>(10, 1));

  auto dense = util::geo::densify(a, 1);

  EXPECT(dense.size() == (size_t)10);

  for (int i = 0; i < 10; i++) {
    EXPECT(dense[i].getX() == approx(i + 1.0));
  }

  dense = util::geo::simplify(dense, 0.1);
  EXPECT(dense.size() == (size_t)2);

  Line<double> b;
  b.push_back(Point<double>(1, 1));
  b.push_back(Point<double>(5, 7));
  b.push_back(Point<double>(10, 3));

  dense = util::geo::densify(b, 1);

  dense = util::geo::simplify(dense, 0.1);
  EXPECT(dense.size() == (size_t)3);
}},

// ___________________________________________________________________________
{
CASE("summed frechet distance") {
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
  EXPECT(fd == approx(2));
}},

// ___________________________________________________________________________
{
CASE("frechet distance") {
  Line<double> e;
  e.push_back(Point<double>(1, 1));
  e.push_back(Point<double>(1, 2));

  Line<double> f;
  f.push_back(Point<double>(1, 1));
  f.push_back(Point<double>(1, 2));

  double fd = util::geo::frechetDist(e, f, 0.1);

  EXPECT(fd == approx(0));

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

  EXPECT(fd == approx(1));

  Line<double> c;
  c.push_back(Point<double>(1, 1));
  c.push_back(Point<double>(2, 1));

  Line<double> d;
  d.push_back(Point<double>(3, 1));
  d.push_back(Point<double>(4, 1));

  fd = util::geo::frechetDist(c, d, 0.1);

  EXPECT(fd == approx(2));


  Line<double> g;
  g.push_back(Point<double>(1, 1));
  g.push_back(Point<double>(10, 1));

  Line<double> h;
  h.push_back(Point<double>(1, 1));
  h.push_back(Point<double>(3, 2));
  h.push_back(Point<double>(3, 1));
  h.push_back(Point<double>(10, 1));

  fd = util::geo::frechetDist(g, h, 0.1);

  EXPECT(fd == approx(1));
}},

// ___________________________________________________________________________
{
CASE("geo box alignment") {
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

  EXPECT(parallelity(box, ml) == approx(1));
  ml = rotate(ml, 45);
  EXPECT(parallelity(box, ml) == approx(0));
  ml = rotate(ml, 45);
  EXPECT(parallelity(box, ml) == approx(1));
  ml = rotate(ml, 45);
  EXPECT(parallelity(box, ml) == approx(0));
  ml = rotate(ml, 45);
  EXPECT(parallelity(box, ml) == approx(1));
}},

// ___________________________________________________________________________
{
CASE("url decode") {
  EXPECT("zürich" == util::urlDecode("z%C3%BCrich"));
  EXPECT("!@$%^*()" == util::urlDecode("!%40%24%25%5E*()"));
  EXPECT("Løkken" == util::urlDecode("L%C3%B8kken"));
  EXPECT("á é" == util::urlDecode("%C3%A1%20%C3%A9"));
  EXPECT("á é" == util::urlDecode("%C3%A1+%C3%A9"));
}},

// ___________________________________________________________________________
{
CASE("json escape") {
  EXPECT("Hello\\\\Goodbye!" == util::jsonStringEscape("Hello\\Goodbye!"));
  EXPECT("\\\"Hello\\\"" == util::jsonStringEscape("\"Hello\""));
}},

// ___________________________________________________________________________
{
CASE("split") {
  EXPECT(util::split("hello,again", ',').size() == (size_t)2);
  EXPECT(util::split("hello,,again", ',').size() == (size_t)3);
  EXPECT(util::split("hello", ',').size() == (size_t)1);
  EXPECT(util::split("", ',').size() == (size_t)0);
}},

// ___________________________________________________________________________
{
CASE("editdist") {
  EXPECT(util::editDist("hello", "mello") == (size_t)1);
  EXPECT(util::editDist("mello", "hello") == (size_t)1);
  EXPECT(util::editDist("abcde", "abfde") == (size_t)1);
  EXPECT(util::editDist("abcd", "abcde") == (size_t)1);
  EXPECT(util::editDist("xabcd", "abcde") == (size_t)2);
  EXPECT(util::editDist("abcd", "abcdes") == (size_t)2);
  EXPECT(util::editDist("hello", "hello") == (size_t)0);
}},

// ___________________________________________________________________________
{
CASE("toString") {
  EXPECT(util::toString(34) == "34");
  EXPECT(util::toString("34") == "34");
}},

// ___________________________________________________________________________
{
CASE("replace") {
  std::string a("lorem ipsum ipsum lorem");

  EXPECT(util::replace(a, "ips", "aa"));
  EXPECT(a == "lorem aaum ipsum lorem");

  EXPECT(!util::replace(a, "blablabla", ""));
  EXPECT(a == "lorem aaum ipsum lorem");

  EXPECT(util::replace(a, "m", ""));
  EXPECT(a == "lore aaum ipsum lorem");

  EXPECT(!util::replace(a, "", ""));
  EXPECT(a == "lore aaum ipsum lorem");

  std::string b("lorem ipsum ipsum lorem");
  EXPECT(util::replaceAll(b, "ips", "aa"));
  EXPECT(b == "lorem aaum aaum lorem");

  EXPECT(util::replaceAll(b, "m", ""));
  EXPECT(b == "lore aau aau lore");

  EXPECT(util::replaceAll(b, "a", "aa"));
  EXPECT(b == "lore aaaau aaaau lore");

  EXPECT(util::replaceAll(b, "e", "e"));
  EXPECT(b == "lore aaaau aaaau lore");

  EXPECT(util::replaceAll(b, "e", "ee"));
  EXPECT(b == "loree aaaau aaaau loree");

  EXPECT(!util::replaceAll(b, "", "ee"));
  EXPECT(b == "loree aaaau aaaau loree");
}},

// ___________________________________________________________________________
{
CASE("Edge-based Dijkstra directed, 1 to all") {
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
    EXPECT(single == u.second);
  }

  // all to 1
  auto eBC = g.addEdg(b, c, 10);

  auto costb = EDijkstra::shortestPathRev(eBC, cFunc);
  for (auto u : costb) {
    int single = EDijkstra::shortestPath(u.first, eBC, cFunc);
    EXPECT(single == u.second);
  }
}},

// ___________________________________________________________________________
{
CASE("Edge-based Dijkstra undirected, edge 1 to 1") {
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

  EXPECT(cost == 2);

  EXPECT(resE.size() == (size_t)3);
  EXPECT(res.size() == (size_t)3);
  EXPECT((*(res.rbegin()))->pl() == "A");
  EXPECT((*(++res.rbegin()))->pl() == "C");
  EXPECT((*(++++res.rbegin()))->pl() == "D");

  EXPECT((*(resE.rbegin())) == eAB);
  EXPECT((*(++resE.rbegin())) == eAC);
  EXPECT((*(++++resE.rbegin())) == eDC);

  cost = EDijkstra::shortestPath(eAB, b, cFunc, &resE, &res);
  EXPECT(cost == 0);
}},

// ___________________________________________________________________________
{
CASE("Edge-based Dijkstra undirected, edge 1 to n") {
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
    int operator()(const Edge<std::string, int>* from, const Node<std::string, int>* n,
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
  EXPECT(cost == 0);
}},

// ___________________________________________________________________________
{
CASE("Edge-based Dijkstra undirected, 1 to n") {
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
    int operator()(const Edge<std::string, int>* from, const Node<std::string, int>* n,
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

  std::unordered_map<Edge<std::string, int>*, EDijkstra::EList<std::string, int>*> resE;
  resE[eDC] = new EDijkstra::EList<std::string, int>();
  resE[eED] = new EDijkstra::EList<std::string, int>();
  std::unordered_map<Edge<std::string, int>*, EDijkstra::NList<std::string, int>*> res;
  res[eDC] = new EDijkstra::NList<std::string, int>();
  res[eED] = new EDijkstra::NList<std::string, int>();
  auto hFunc = EDijkstra::ZeroHeurFunc<std::string, int, int>();
  std::unordered_map<Edge<std::string, int>*, int> cost = EDijkstra::shortestPath(eAB, tos, cFunc, hFunc, resE, res);

  EXPECT(cost[eDC] == 2);
  EXPECT(cost[eED] == 2);

  EXPECT(resE[eDC]->size() == (size_t)3);
  EXPECT(res[eED]->size() == (size_t)3);

  EXPECT(resE[eDC]->size() == (size_t)3);
  EXPECT(res[eED]->size() == (size_t)3);
}},

// ___________________________________________________________________________
{
CASE("Edge-based Dijkstra undirected") {
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
    int operator()(const Edge<std::string, int>* fr, const Node<std::string, int>* n,
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

  EXPECT(res.size() == (size_t)5);
  EXPECT((*(res.rbegin()))->pl() == "A");
  EXPECT((*(++res.rbegin()))->pl() == "C");
  EXPECT((*(++++res.rbegin()))->pl() == "D");
  EXPECT((*(++++++res.rbegin()))->pl() == "E");
  EXPECT((*(++++++++res.rbegin()))->pl() == "B");
  EXPECT(cost == 4);
  EXPECT((*(resE.rbegin()))->getFrom()->pl() == "A");
  EXPECT((*(++resE.rbegin()))->getFrom()->pl() == "D");
  EXPECT((*(++++resE.rbegin()))->getFrom()->pl() == "E");
  EXPECT((*(++++++resE.rbegin()))->getTo()->pl() == "B");

  EXPECT(resE.size() == (size_t)4);

  cost = EDijkstra::shortestPath(d, b, cFunc, &res);
  EXPECT(cost == 2);

  cost = EDijkstra::shortestPath(b, d, cFunc, &res);
  EXPECT(cost == 2);

  cost = EDijkstra::shortestPath(e, b, cFunc, &res);
  EXPECT(cost == 1);

  cost = EDijkstra::shortestPath(b, e, cFunc, &res);
  EXPECT(cost == 1);

  cost = EDijkstra::shortestPath(b, a, cFunc, &res);
  EXPECT(cost == 4);

  cost = EDijkstra::shortestPath(c, a, cFunc, &res);
  EXPECT(cost == 1);

  cost = EDijkstra::shortestPath(a, c, cFunc, &res);
  EXPECT(cost == 1);

  cost = EDijkstra::shortestPath(a, d, cFunc, &res);
  EXPECT(cost == 2);
}},

// ___________________________________________________________________________
{
CASE("Edge-based Dijkstra") {
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

  EXPECT(cost == 3);

  g.addEdg(c, d, 3);
  cost = EDijkstra::shortestPath(a, d, cFunc, &res);

  EXPECT(cost == 3);

  g.addEdg(a, b, 1);
  g.addEdg(x, a, 1);
  cost = EDijkstra::shortestPath(a, d, cFunc, &res);

  EXPECT(cost == 2);
}},

// ___________________________________________________________________________
{
CASE("Dijkstra") {
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

  EXPECT(cost == 3);
  EXPECT(res.size() == (size_t)4);

  g.addEdg(c, d, 3);
  cost = Dijkstra::shortestPath(a, d, cFunc, &res);

  EXPECT(cost == 3);

  g.addEdg(a, b, 1);
  g.addEdg(x, a, 1);
  cost = Dijkstra::shortestPath(a, d, cFunc, &res);

  EXPECT(cost == 2);

  const std::set<Node<int, int>*> to{b, c, d, x};
  std::unordered_map<Node<int, int>*, Dijkstra::EList<int, int>*> resEdges;
  std::unordered_map<Node<int, int>*, Dijkstra::NList<int, int>*> resNodes;

  for (auto n : to) {
    resEdges[n] = new Dijkstra::EList<int, int>();
    resNodes[n] = new Dijkstra::NList<int, int>();
  }

  auto costs = Dijkstra::shortestPath(a, to, cFunc, resEdges, resNodes);

  EXPECT(costs[b] == 1);
  EXPECT(costs[c] == 1);
  EXPECT(costs[d] == 2);
  EXPECT(costs[x] == 999);
}},

// ___________________________________________________________________________
{
CASE("nullable") {
  {
    util::Nullable<std::string> nullable;
    EXPECT(nullable.isNull());
  }

  {
    util::Nullable<std::string> nullable(0);
    EXPECT(nullable.isNull());
  }

  {
    std::string str = "aa";
    util::Nullable<std::string> nullable(&str);
    EXPECT(!nullable.isNull());

    EXPECT(nullable == "aa");
    EXPECT(!(nullable == "aaa"));
    EXPECT(!(nullable != "aa"));
    EXPECT(nullable == "aa");

    EXPECT(nullable.get() == "aa");
    EXPECT(std::string(nullable) == "aa");
  }

  {
    int a = 23;
    util::Nullable<int> nullable(a);
    util::Nullable<int> nullable2(24);
    EXPECT(!nullable.isNull());

    EXPECT(nullable == 23);
    EXPECT(nullable >= 23);
    EXPECT(nullable <= 23);
    EXPECT(nullable < 24);
    EXPECT(nullable < 24);
    EXPECT(!(nullable < 22));
    EXPECT(nullable != nullable2);
    EXPECT(nullable < nullable2);
    EXPECT(nullable2 > nullable);

    util::Nullable<int> nullable3(nullable);
    EXPECT(nullable == nullable3);

    nullable3 = nullable2;
    EXPECT(nullable2 == nullable3);
    EXPECT(nullable3 == 24);
    EXPECT(nullable2 == 24);
    EXPECT(nullable2 == nullable2.get());
    EXPECT(int(nullable2) == nullable2.get());
    EXPECT(!nullable3.isNull());
    EXPECT(!nullable2.isNull());

    util::Nullable<int> voidnull;
    EXPECT(voidnull.isNull());

    EXPECT_THROWS(nullable == voidnull);
  }
}},

// ___________________________________________________________________________
{
CASE("geometry") {
  geo::Point<double> a(1, 2);
  geo::Point<double> b(2, 3);
  geo::Point<double> c(4, 5);
  EXPECT(a.getX() == approx(1));
  EXPECT(a.getY() == approx(2));

  a.setX(3);
  EXPECT(a.getX() == approx(3));
  EXPECT(a.getY() == approx(2));

  a.setY(4);
  EXPECT(a.getX() == approx(3));
  EXPECT(a.getY() == approx(4));

  auto d = a + b;
  EXPECT(d.getX() == approx(5));
  EXPECT(d.getY() == approx(7));

  a.setX(1);
  a.setY(2);

  EXPECT(geo::dist(a, a) == approx(0));
  EXPECT(geo::dist(a, b) == approx(sqrt(2)));

  d = d + d;

  geo::Box<double> box(a, c);
  EXPECT(geo::contains(a, box));
  EXPECT(geo::contains(b, box));
  EXPECT(geo::contains(c, box));
  EXPECT(!geo::contains(d, box));

  geo::Line<double> line{a, b, c};

  EXPECT(geo::contains(line, box));
  line.push_back(d);
  EXPECT(!geo::contains(line, box));

  geo::LineSegment<double> ls{a, b};
  EXPECT(geo::contains(a, ls));
  EXPECT(geo::contains(b, ls));
  EXPECT(!geo::contains(c, ls));
  EXPECT(geo::contains(a + geo::Point<double>(.5, .5), ls));
  EXPECT(!geo::contains(a + geo::Point<double>(1.5, 1.5), ls));

  geo::LineSegment<double> lsa{geo::Point<double>(1, 1), geo::Point<double>(2, 2)};
  geo::LineSegment<double> lsb{geo::Point<double>(1, 2), geo::Point<double>(2, 1)};
  geo::LineSegment<double> lsc{geo::Point<double>(2.1, 2), geo::Point<double>(3, 3)};

  EXPECT(geo::crossProd(lsa.first, lsb) == approx(-1));
  EXPECT(geo::crossProd(lsa.second, lsb) == approx(1));

  EXPECT(geo::intersects(lsa, lsb));

  EXPECT(geo::intersects(lsa, lsa));
  EXPECT(geo::intersects(lsb, lsb));
  EXPECT(!geo::intersects(lsa, lsc));

  EXPECT(!geo::intersects(geo::Point<double>(871569.2, 6104550.4), geo::Point<double>(871581.2, 6104536), geo::Point<double>(871580.3, 6104541.3), geo::Point<double>(871625.7, 6104510.1)));

  geo::Line<double> l{geo::Point<double>(1, 1), geo::Point<double>(2, 2), geo::Point<double>(2, 4)};
  EXPECT(!geo::contains(geo::Point<double>(1, 2), l));
  EXPECT(geo::contains(geo::Point<double>(2, 2), l));
  EXPECT(geo::contains(geo::Point<double>(2, 3), l));

  geo::Box<double> bbox(geo::Point<double>(1, 1), geo::Point<double>(3, 3));
  EXPECT(geo::intersects(l, bbox));
  geo::Line<double> ll{geo::Point<double>(0, 0), geo::Point<double>(4, 4)};
  EXPECT(geo::intersects(ll, bbox));
  geo::Line<double> lll{geo::Point<double>(0, 0), geo::Point<double>(0, 4)};
  EXPECT(!geo::intersects(lll, bbox));
  geo::Line<double> llll{geo::Point<double>(1.2, 0), geo::Point<double>(1, 2)};
  EXPECT(geo::intersects(llll, bbox));

  Line<double> l5new;
  l5new.push_back(Point<double>(-10, -5));
  l5new.push_back(Point<double>(-8, -4));
  EXPECT(geo::getBoundingBox(l5new).getUpperRight().getX() == approx(-8));
  EXPECT(geo::getBoundingBox(l5new).getUpperRight().getY() == approx(-4));

  Line<double> l5;
  l5.push_back(Point<double>(0, 0));
  l5.push_back(Point<double>(1.5, 2));
  Box<double> req(Point<double>(.5, 1), Point<double>(1, 1.5));

  EXPECT(geo::getBoundingBox(l5[0]).getLowerLeft().getX() == approx(0));
  EXPECT(geo::getBoundingBox(l5[0]).getLowerLeft().getY() == approx(0));

  EXPECT(geo::getBoundingBox(l5).getLowerLeft().getX() == approx(0));
  EXPECT(geo::getBoundingBox(l5).getLowerLeft().getY() == approx(0));
  EXPECT(geo::getBoundingBox(l5).getUpperRight().getX() == approx(1.5));
  EXPECT(geo::getBoundingBox(l5).getUpperRight().getY() == approx(2));
  EXPECT(geo::intersects(geo::getBoundingBox(l5), geo::getBoundingBox(Line<double>{Point<double>(.5, 1), Point<double>(1, 1)})));
  EXPECT(geo::intersects(l5, Line<double>{Point<double>(.5, 1), Point<double>(1, 1)}));
  EXPECT(geo::intersects(l5, req));

  Box<double> boxa(Point<double>(1, 1), Point<double>(2, 2));
  EXPECT(geo::intersects(boxa, Box<double>(Point<double>(1.5, 1.5), Point<double>(1.7, 1.7))));
  EXPECT(geo::intersects(boxa, Box<double>(Point<double>(0, 0), Point<double>(3, 3))));
  EXPECT(geo::intersects(boxa, Box<double>(Point<double>(1.5, 1.5), Point<double>(3, 3))));
  EXPECT(geo::intersects(boxa, Box<double>(Point<double>(0, 0), Point<double>(1.5, 1.5))));

  EXPECT(geo::intersects(Box<double>(Point<double>(1.5, 1.5), Point<double>(1.7, 1.7)), boxa));
  EXPECT(geo::intersects(Box<double>(Point<double>(0, 0), Point<double>(3, 3)), boxa));
  EXPECT(geo::intersects(Box<double>(Point<double>(1.5, 1.5), Point<double>(3, 3)), boxa));
  EXPECT(geo::intersects(Box<double>(Point<double>(0, 0), Point<double>(1.5, 1.5)), boxa));

  Polygon<double> poly({Point<double>(1, 1), Point<double>(3, 2), Point<double>(4, 3), Point<double>(6, 3), Point<double>(5, 1)});
  EXPECT(geo::getWKT(poly) == "POLYGON ((1 1, 3 2, 4 3, 6 3, 5 1, 1 1))");
  EXPECT(geo::contains(Point<double>(4, 2), poly));
  EXPECT(!geo::contains(Point<double>(3, 3), poly));
  EXPECT(geo::contains(Point<double>(1, 1), poly));
  EXPECT(geo::contains(Point<double>(3, 2), poly));
  EXPECT(geo::contains(Point<double>(4, 3), poly));
  EXPECT(geo::contains(Point<double>(6, 3), poly));
  EXPECT(geo::contains(Point<double>(5, 1), poly));

  EXPECT(geo::contains(Line<double>{Point<double>(6, 3), Point<double>(5, 1)}, poly));
  EXPECT(!geo::contains(Line<double>{Point<double>(6, 3), Point<double>(50, 1)}, poly));
  EXPECT(geo::contains(Line<double>{Point<double>(4, 2), Point<double>(4.5, 2)}, poly));
  EXPECT(geo::contains(Line<double>{Point<double>(4, 2), Point<double>(5, 1)}, poly));

  Box<double> polybox(Point<double>(1, 1), Point<double>(6, 4));
  EXPECT(geo::centroid(polybox).getX() == approx(3.5));
  EXPECT(geo::centroid(polybox).getY() == approx(2.5));
  EXPECT(geo::contains(poly, polybox));
  EXPECT(!geo::contains(polybox, poly));
  Box<double> polybox2(Point<double>(4, 1), Point<double>(5, 2));
  EXPECT(geo::contains(polybox2, poly));
  EXPECT(geo::contains(poly, getBoundingBox(poly)));

  Point<double> rotP(2, 2);
  EXPECT(geo::dist(geo::rotate(rotP, 180, Point<double>(1, 1)), Point<double>(0, 0)) == approx(0));
  EXPECT(geo::dist(geo::rotate(rotP, 360, Point<double>(1, 1)), rotP) == approx(0));

  Line<double> rotLine({{1, 1}, {3, 3}});
  EXPECT(geo::rotate(rotLine, 90, Point<double>(2, 2))[0].getX() == approx(1));
  EXPECT(geo::rotate(rotLine, 90, Point<double>(2, 2))[0].getY() == approx(3));
  EXPECT(geo::rotate(rotLine, 90, Point<double>(2, 2))[1].getX() == approx(3));
  EXPECT(geo::rotate(rotLine, 90, Point<double>(2, 2))[1].getY() == approx(1));

  MultiLine<double> multiRotLine({{{1, 1}, {3, 3}}, {{1, 3}, {3, 1}}});
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][0].getX() == approx(1));
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][0].getY() == approx(3));
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][1].getX() == approx(3));
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[0][1].getY() == approx(1));
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][0].getX() == approx(3));
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][0].getY() == approx(3));
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][1].getX() == approx(1));
  EXPECT(geo::rotate(multiRotLine, 90, Point<double>(2, 2))[1][1].getY() == approx(1));

  EXPECT(geo::getWKT(multiRotLine) == "MULTILINESTRING ((1 1, 3 3), (1 3, 3 1))");

  EXPECT(geo::contains(multiRotLine[0], geo::move(geo::move(multiRotLine, 1.0, 2.0), -1.0, -2.0)[0]));
  EXPECT(geo::contains(multiRotLine, geo::getBoundingBox(Line<double>{{1, 1}, {3, 3}, {1, 3}, {3, 1}})));

  EXPECT(geo::contains(getBoundingBox(multiRotLine), geo::getBoundingBox(Line<double>{{1, 1}, {3, 3}, {1, 3}, {3, 1}})));
  EXPECT(geo::contains(geo::getBoundingBox(Line<double>{{1, 1}, {3, 3}, {1, 3}, {3, 1}}), getBoundingBox(multiRotLine)));

  EXPECT(geo::dist(geo::centroid(rotP), rotP) == approx(0));
  EXPECT(geo::dist(geo::centroid(rotLine), rotP) == approx(0));
  EXPECT(geo::dist(geo::centroid(polybox), Point<double>(3.5, 2.5)) == approx(0));
  EXPECT(geo::dist(geo::centroid(Polygon<double>({{0, 0}, {3, 4}, {4,3}})), Point<double>(7.0/3.0,7.0/3.0)) == approx(0));

  auto polyy = Polygon<double>({{0, 0}, {3, 4}, {4,3}});
  MultiPolygon<double> mpoly{polyy, polyy};

  EXPECT(geo::getWKT(polyy) == "POLYGON ((0 0, 3 4, 4 3, 0 0))");
  EXPECT(geo::getWKT(mpoly) == "MULTIPOLYGON (((0 0, 3 4, 4 3, 0 0)), ((0 0, 3 4, 4 3, 0 0)))");

  auto hull = geo::convexHull(Line<double>{{0.1, 3}, {1, 1}, {2, 2}, {4, 4}, {0, 0}, {1, 2}, {3, 1}, {3, 3}});
  EXPECT(hull.getOuter().size() == size_t(4));
  EXPECT(hull.getOuter()[0].getX() == approx(0));
  EXPECT(hull.getOuter()[0].getY() == approx(0));
  EXPECT(hull.getOuter()[1].getX() == approx(3));
  EXPECT(hull.getOuter()[1].getY() == approx(1));
  EXPECT(hull.getOuter()[2].getX() == approx(4));
  EXPECT(hull.getOuter()[2].getY() == approx(4));
  EXPECT(hull.getOuter()[3].getX() == approx(0.1));
  EXPECT(hull.getOuter()[3].getY() == approx(3));
  EXPECT(geo::contains(geo::convexHull(geo::getBoundingBox(poly)), geo::getBoundingBox(poly)));
  EXPECT(geo::contains(geo::getBoundingBox(poly), geo::convexHull(geo::getBoundingBox(poly))));

  auto hull2 = geo::convexHull(Line<double>{{0.1, 3}, {1, 1}, {2, 2}, {4, 4}, {0, 0}, {1, 2}, {3, 1}, {3, 3}, {-0.1, 1}});
  EXPECT(hull2.getOuter().size() == size_t(5));
  EXPECT(hull2.getOuter()[0].getX() == approx(-.1));
  EXPECT(hull2.getOuter()[0].getY() == approx(1));
  EXPECT(hull2.getOuter()[1].getX() == approx(0));
  EXPECT(hull2.getOuter()[1].getY() == approx(0));
  EXPECT(hull2.getOuter()[2].getX() == approx(3));
  EXPECT(hull2.getOuter()[2].getY() == approx(1));
  EXPECT(hull2.getOuter()[3].getX() == approx(4));
  EXPECT(hull2.getOuter()[3].getY() == approx(4));
  EXPECT(hull2.getOuter()[4].getX() == approx(0.1));
  EXPECT(hull2.getOuter()[4].getY() == approx(3));

  auto hull3 = geo::convexHull(Line<double>{{0.1, 3}, {4, 4}, {0, 0}, {1, 2}, {3, 1}});
  EXPECT(hull3.getOuter().size() == size_t(4));
  EXPECT(hull3.getOuter()[0].getX() == approx(0));
  EXPECT(hull3.getOuter()[0].getY() == approx(0));
  EXPECT(hull3.getOuter()[1].getX() == approx(3));
  EXPECT(hull3.getOuter()[1].getY() == approx(1));
  EXPECT(hull3.getOuter()[2].getX() == approx(4));
  EXPECT(hull3.getOuter()[2].getY() == approx(4));
  EXPECT(hull3.getOuter()[3].getX() == approx(0.1));
  EXPECT(hull3.getOuter()[3].getY() == approx(3));

  hull3 = geo::convexHull(Line<double>{{0.1, 3}, {4, 4}, {2, 1}, {3, 2}, {0, 0}, {1, 2}, {3, 1}});
  EXPECT(hull3.getOuter().size() == size_t(4));
  EXPECT(hull3.getOuter()[0].getX() == approx(0));
  EXPECT(hull3.getOuter()[0].getY() == approx(0));
  EXPECT(hull3.getOuter()[1].getX() == approx(3));
  EXPECT(hull3.getOuter()[1].getY() == approx(1));
  EXPECT(hull3.getOuter()[2].getX() == approx(4));
  EXPECT(hull3.getOuter()[2].getY() == approx(4));
  EXPECT(hull3.getOuter()[3].getX() == approx(0.1));
  EXPECT(hull3.getOuter()[3].getY() == approx(3));

  hull3 = geo::convexHull(Line<double>{{4, 4}, {1, 2}, {2, 1}, {3, 2}, {0.1, 3}, {0, 0}, {1, 2}, {3, 1}});
  EXPECT(hull3.getOuter().size() == size_t(4));
  EXPECT(hull3.getOuter()[0].getX() == approx(0));
  EXPECT(hull3.getOuter()[0].getY() == approx(0));
  EXPECT(hull3.getOuter()[1].getX() == approx(3));
  EXPECT(hull3.getOuter()[1].getY() == approx(1));
  EXPECT(hull3.getOuter()[2].getX() == approx(4));
  EXPECT(hull3.getOuter()[2].getY() == approx(4));
  EXPECT(hull3.getOuter()[3].getX() == approx(0.1));
  EXPECT(hull3.getOuter()[3].getY() == approx(3));

  hull3 = geo::convexHull(Line<double>{{4, 4}, {1, 2}, {3, 1}});
  EXPECT(hull3.getOuter().size() == size_t(3));
  EXPECT(hull3.getOuter()[0].getX() == approx(1));
  EXPECT(hull3.getOuter()[0].getY() == approx(2));
  EXPECT(hull3.getOuter()[1].getX() == approx(3));
  EXPECT(hull3.getOuter()[1].getY() == approx(1));
  EXPECT(hull3.getOuter()[2].getX() == approx(4));
  EXPECT(hull3.getOuter()[2].getY() == approx(4));

  hull3 = geo::convexHull(Line<double>{{4, 4}, {1, 2}, {3, 10}});
  EXPECT(hull3.getOuter().size() == size_t(3));
  EXPECT(hull3.getOuter()[0].getX() == approx(1));
  EXPECT(hull3.getOuter()[0].getY() == approx(2));
  EXPECT(hull3.getOuter()[1].getX() == approx(4));
  EXPECT(hull3.getOuter()[1].getY() == approx(4));
  EXPECT(hull3.getOuter()[2].getX() == approx(3));
  EXPECT(hull3.getOuter()[2].getY() == approx(10));

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
  EXPECT(geo::contains(test, hull3));
  EXPECT(hull3.getOuter().size() == size_t(8));
  EXPECT(geo::contains(Polygon<double>({{-0.161920957418085, -0.4055339716426413},
                                        {0.05054295812784038, 0.4754929463150845},
                                        {0.4823896228171788, -0.4776170002088109},
                                        {0.4932166845474547, 0.4928094162538735},
                                        {-0.3521487911717489, 0.4352656197131292},
                                        {-0.4907368011686362, 0.1865826865533206},
                                        {0.4916198379282093, -0.345391701297268},
                                        {-0.4404289572876217,
                                         -0.2894855991839297}}), hull3));
  EXPECT(geo::contains(hull3, Polygon<double>({{-0.161920957418085, -0.4055339716426413},
                                        {0.05054295812784038, 0.4754929463150845},
                                        {0.4823896228171788, -0.4776170002088109},
                                        {0.4932166845474547, 0.4928094162538735},
                                        {-0.3521487911717489, 0.4352656197131292},
                                        {-0.4907368011686362, 0.1865826865533206},
                                        {0.4916198379282093, -0.345391701297268},
                                        {-0.4404289572876217,
                                         -0.2894855991839297}})));

  hull3 = geo::convexHull(Line<double>{{3, 6}, {8, 10}, {3, 5}, {20, -10}, {-4, 5}, {10, 2}, {5, 1}, {45, 1}, {30, -9}, {3, 14}, {25, -5.5}});
  EXPECT(hull3.getOuter().size() == size_t(5));
  EXPECT(hull3.getOuter()[0].getX() == approx(-4));
  EXPECT(hull3.getOuter()[0].getY() == approx(5));
  EXPECT(hull3.getOuter()[1].getX() == approx(20));
  EXPECT(hull3.getOuter()[1].getY() == approx(-10));
  EXPECT(hull3.getOuter()[2].getX() == approx(30));
  EXPECT(hull3.getOuter()[2].getY() == approx(-9));
  EXPECT(hull3.getOuter()[3].getX() == approx(45));
  EXPECT(hull3.getOuter()[3].getY() == approx(1));
  EXPECT(hull3.getOuter()[4].getX() == approx(3));
  EXPECT(hull3.getOuter()[4].getY() == approx(14));

  hull3 = geo::convexHull(Line<double>{{7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}});
  EXPECT(hull3.getOuter().size() == size_t(8));
  EXPECT(geo::contains(geo::Polygon<double>({{-9, 0}, {-7, -7}, {0, -9}, {7, -7}, {9, 0}, {7, 7}, {0, 9}, {-7, 7}}), hull3));
  EXPECT(geo::contains(hull3, geo::Polygon<double>({{-9, 0}, {-7, -7}, {0, -9}, {7, -7}, {9, 0}, {7, 7}, {0, 9}, {-7, 7}})));

  hull3 = geo::convexHull(Line<double>{{7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}, {0, 0}, {1, 2}, {-2, 1}, {-1, -1}, {3, 4}, {4, 3}, {-5, 4}, {6, 5}});
  EXPECT(hull3.getOuter().size() == size_t(8));
  EXPECT(geo::contains(geo::Polygon<double>({{-9, 0}, {-7, -7}, {0, -9}, {7, -7}, {9, 0}, {7, 7}, {0, 9}, {-7, 7}}), hull3));
  EXPECT(geo::contains(hull3, geo::Polygon<double>({{-9, 0}, {-7, -7}, {0, -9}, {7, -7}, {9, 0}, {7, 7}, {0, 9}, {-7, 7}})));

  hull3 = geo::convexHull(Line<double>{{0, 0}, {1, 2}, {-2, 1}, {-1, -1}, {3, 4}, {4, 3}, {-5, 4}, {6, 5}, {7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}, {-8, 0}, {8, 0}, {-7, 0}, {7, 0}, {-6, 0}, {6, 0}, {-5, 0}, {5, 0}, {-4, 0}, {4, 0}, {-3, 0}, {3, 0}, {-2, 0}, {2, 0}, {-1, 0}, {1, 0}, {0, -8}, {0, 8}, {0, -7}, {0, 7}, {0, -6}, {0, 6}, {0, -5}, {0, 5}, {0, -4}, {0, 4}, {0, -3}, {0, 3}, {0, -2}, {0, 2}, {0, -1}, {0, 1}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}, {1, -1}, {2, -2}, {3, -3}, {4, -4},  {5, -5}, {6, -6}, {-1, 1}, {-2, 2}, {-3, 3}, {-4, 4}, {-5, 5}, {-6, 6}, {-1, -1}, {-2, -2}, {-3, -3}, {-4, -4}, {-5, -5}, {-6, -6}});
  EXPECT(hull3.getOuter().size() == size_t(8));
  EXPECT(geo::contains(geo::Polygon<double>({{-9, 0}, {-7, -7}, {0, -9}, {7, -7}, {9, 0}, {7, 7}, {0, 9}, {-7, 7}}), hull3));
  EXPECT(geo::contains(hull3, geo::Polygon<double>({{-9, 0}, {-7, -7}, {0, -9}, {7, -7}, {9, 0}, {7, 7}, {0, 9}, {-7, 7}})));

  EXPECT(geo::area(geo::Point<double>(1, 2)) == approx(0));
  EXPECT(geo::area(geo::Line<double>{{1, 2}, {2, 5}}) == approx(0));
  EXPECT(geo::area(geo::Box<double>({0, 0}, {1, 1})) == approx(1));
  EXPECT(geo::area(geo::Box<double>({1, 1}, {1, 1})) == approx(0));
  EXPECT(geo::area(geo::Box<double>({0, 0}, {2, 2})) == approx(4));
  EXPECT(geo::area(geo::Polygon<double>({{0, 0}, {1, 0}, {1, 1}, {0, 1}})) == approx(1));
  EXPECT(geo::area(geo::Polygon<double>({{0, 0}, {1, 0}, {1, 1}})) == approx(0.5));

  auto obox = geo::getOrientedEnvelope(geo::Line<double>{{0, 0}, {1, 1}, {1.5, 0.5}});
  EXPECT(geo::contains(geo::convexHull(obox), geo::Polygon<double>({{0.0, 0.0}, {1.0, 1.0}, {1.5, 0.5}, {0.5, -0.5}})));
  EXPECT(geo::contains(geo::Polygon<double>({{0.0, 0.0}, {1.0, 1.0}, {1.5, 0.5}, {0.5, -0.5}}), geo::convexHull(obox)));

  EXPECT(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}}, geo::LineSegment<double>{{2, 2}, {2, 0}}) == approx(0));
  EXPECT(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}}, geo::LineSegment<double>{{2, 4}, {2, 2}}) == approx(1));
  EXPECT(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}}, geo::LineSegment<double>{{1, 1}, {3, 1}}) == approx(0));
  EXPECT(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}}, geo::LineSegment<double>{{1, 2}, {3, 2}}) == approx(1));
  EXPECT(geo::dist(geo::LineSegment<double>{{1, 1}, {3, 1}}, geo::LineSegment<double>{{1, 2}, {3, 5}}) == approx(1));

  EXPECT(geo::dist(geo::Line<double>{{1, 1}, {3, 1}}, geo::Point<double>{2, 1}) == approx(0));
  EXPECT(geo::dist(geo::Line<double>{{1, 1}, {3, 1}}, geo::Point<double>{2, 2}) == approx(1));
  EXPECT(geo::dist(geo::Line<double>{{1, 1}, {3, 1}}, geo::Point<double>{3, 1}) == approx(0));
  EXPECT(geo::dist(geo::Line<double>{{1, 1}, {3, 1}}, geo::Point<double>{1, 1}) == approx(0));

  EXPECT(geo::dist(Line<double>{{7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}}, Line<double>{{7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}}) == approx(0));
  EXPECT(geo::dist(Line<double>{{7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}}, LineSegment<double>{{6, 7}, {8, -7}}) == approx(0));
  EXPECT(geo::dist(Line<double>{{7, 7}, {7, -7}, {-7, -7}, {-7, 7}, {9, 0}, {-9, 0}, {0, 9}, {0, -9}}, Point<double>{7, 4}) == approx(0));
  EXPECT(geo::dist(Line<double>{{0, 0}, {1, 1}, {2, 0}}, Line<double>{{1.5, 0.5}, {1.5, 100}}) == approx(0));
  EXPECT(geo::dist(Line<double>{{0, 0}, {1, 1}, {2, 0}}, Line<double>{{2, 0.5}, {2, 100}}) == approx(0.353553));
}

}};

// _____________________________________________________________________________
int main(int argc, char** argv) {
  return(lest::run(specification, argc, argv));
}
