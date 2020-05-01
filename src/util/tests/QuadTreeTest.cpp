// Copyright 2016
// Author: Patrick Brosi

#include "util/Misc.h"
#include "util/geo/QuadTree.h"
#include "util/tests/QuadTreeTest.h"

using util::approx;
using util::geo::QuadTree;
using util::geo::DPoint;
using util::geo::DBox;

// _____________________________________________________________________________
void QuadTreeTest::run() {
  // ___________________________________________________________________________
  {
    QuadTree<int, double> qt(4, 4, DBox(DPoint(0, 0), DPoint(10, 10)));

    qt.insert(0, {2, 2});
    TEST(qt.size(), ==, 1);

    qt.insert(666, {-1, 0});
    TEST(qt.size(), ==, 1);

    qt.insert(1, {0, 0});
    TEST(qt.size(), ==, 2);

    qt.insert(2, {0, 1});
    TEST(qt.size(), ==, 3);

    qt.insert(3, {6, 9});
    TEST(qt.size(), ==, 4);

    qt.insert(4, {9, 0});
    TEST(qt.size(), ==, 5);
  }

}
