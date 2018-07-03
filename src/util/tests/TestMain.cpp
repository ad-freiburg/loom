// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include "lest.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/Geo.h"
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
CASE("atof") {
  EXPECT(util::atof("45.534215") == approx(45.534215));
  EXPECT(util::atof("5.534") == approx(5.534));
  EXPECT(util::atof("534") == approx(534));
  EXPECT(util::atof("-534") == approx(-534));
  EXPECT(util::atof("-45.534215") == approx(-45.534215));
  EXPECT(util::atof("-45.534215", 2) == approx(-45.53));


  // TODO: more test cases
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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

},

// ___________________________________________________________________________
CASE("grid") {
  Grid<int, Line, double> g(.5, .5, Box<double>(Point<double>(0, 0), Point<double>(3, 3)));

  Line<double> l;
  l.push_back(Point<double>(0, 0));
  l.push_back(Point<double>(1.5, 2));

  Line<double> l2;
  l.push_back(Point<double>(2.5, 1));
  l.push_back(Point<double>(2.5, 2));

  g.add(l, 1);
  g.add(l2, 2);

  std::set<int> ret;

  Box<double> req(Point<double>(.5, 1), Point<double>(1, 1.5));
  g.get(req, &ret);
  EXPECT(ret.size() == (size_t)1);

  ret.clear();
  g.getNeighbors(1, 0, &ret);
  // TODO!
  //EXPECT(ret.size() == 1);

  ret.clear();
  g.getNeighbors(1, 0.55, &ret);
  EXPECT(ret.size() == (size_t)2);


  // TODO: more test cases
},

// ___________________________________________________________________________
CASE("densify") {
  Line<double> a;
  a.push_back(Point<double>(1, 1));
  a.push_back(Point<double>(10, 1));

  auto dense = util::geo::densify(a, 1);

  EXPECT(dense.size() == (size_t)10);

  for (int i = 0; i < 10; i++) {
    EXPECT(dense[i].get<0>() == approx(i + 1.0));
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
CASE("url decode") {
  EXPECT("zürich" == util::urlDecode("z%C3%BCrich"));
  EXPECT("!@$%^*()" == util::urlDecode("!%40%24%25%5E*()"));
  EXPECT("Løkken" == util::urlDecode("L%C3%B8kken"));
  EXPECT("á é" == util::urlDecode("%C3%A1%20%C3%A9"));
  EXPECT("á é" == util::urlDecode("%C3%A1+%C3%A9"));
},

// ___________________________________________________________________________
CASE("json escape") {
  EXPECT("Hello\\\\Goodbye!" == util::jsonStringEscape("Hello\\Goodbye!"));
  EXPECT("\\\"Hello\\\"" == util::jsonStringEscape("\"Hello\""));
},

// ___________________________________________________________________________
CASE("split") {
  EXPECT(util::split("hello,again", ',').size() == (size_t)2);
  EXPECT(util::split("hello,,again", ',').size() == (size_t)3);
  EXPECT(util::split("hello", ',').size() == (size_t)1);
  EXPECT(util::split("", ',').size() == (size_t)0);
},

// ___________________________________________________________________________
CASE("editdist") {
  EXPECT(util::editDist("hello", "mello") == (size_t)1);
  EXPECT(util::editDist("mello", "hello") == (size_t)1);
  EXPECT(util::editDist("abcde", "abfde") == (size_t)1);
  EXPECT(util::editDist("abcd", "abcde") == (size_t)1);
  EXPECT(util::editDist("xabcd", "abcde") == (size_t)2);
  EXPECT(util::editDist("abcd", "abcdes") == (size_t)2);
  EXPECT(util::editDist("hello", "hello") == (size_t)0);
},

// ___________________________________________________________________________
CASE("toString") {
  EXPECT(util::toString(34) == "34");
  EXPECT(util::toString("34") == "34");
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
},

// ___________________________________________________________________________
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
}

};

// _____________________________________________________________________________
int main(int argc, char** argv) {
  return(lest::run(specification, argc, argv));
}
