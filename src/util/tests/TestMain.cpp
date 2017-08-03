// Copyright 2016
// Author: Patrick Brosi
//

#include <string>
#include "lest.h"
#include "util/Nullable.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/graph/Graph.h"
#include "util/geo/Grid.h"

using lest::approx;
using namespace util::geo;
using namespace util::graph;

// define LEST cases
const lest::test specification[] = {

// ___________________________________________________________________________
CASE("graph") {
  Graph<int, int> g;

  Node<int, int>* a = new Node<int, int>(0);
  Node<int, int>* b = new Node<int, int>(0);
  g.addNode(a);
  EXPECT(g.getNodes()->size() == 1);
  g.addNode(b);
  EXPECT(g.getNodes()->size() == 2);


  // TODO: more test cases 
},

// ___________________________________________________________________________
CASE("grid") {
  Grid<int, Line> g(.5, .5, Box(Point(0, 0), Point(3, 3)));

  Line l;
  l.push_back(Point(0, 0));
  l.push_back(Point(1.5, 2));

  Line l2;
  l.push_back(Point(2.5, 1));
  l.push_back(Point(2.5, 2));

  g.add(l, 1);
  g.add(l2, 2);

  std::set<int> ret;

  Box req(Point(.5, 1), Point(1, 1.5));
  g.get(req, &ret);
  EXPECT(ret.size() == 1);

  ret.clear();
  g.getNeighbors(1, 0, &ret);
  EXPECT(ret.size() == 1);

  ret.clear();
  g.getNeighbors(1, 0.55, &ret);
  EXPECT(ret.size() == 2);


  // TODO: more test cases 
},

// ___________________________________________________________________________
CASE("geo box alignment") {
  Line a;
  a.push_back(Point(1, 1));
  a.push_back(Point(1, 2));

  Line b;
  b.push_back(Point(1, 2));
  b.push_back(Point(2, 2));

  Line c;
  c.push_back(Point(2, 2));
  c.push_back(Point(2, 1));

  Line d;
  d.push_back(Point(2, 1));
  d.push_back(Point(1, 1));

  Box box(Point(2, 3), Point(5, 4));
  MultiLine ml;
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
