// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <climits>
#include "shared/linegraph/LineGraph.h"
#include "topo/mapconstructor/MapConstructor.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/log/Log.h"

using topo::MapConstructor;
using topo::ShrdSegWrap;
using topo::config::TopoConfig;

using util::geo::Box;
using util::geo::DBox;
using util::geo::DPoint;
using util::geo::extendBox;
using util::geo::Grid;
using util::geo::Point;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePair;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::Station;
using shared::linegraph::LineOcc;

// _____________________________________________________________________________
MapConstructor::MapConstructor(const TopoConfig* cfg, LineGraph* g)
    : _cfg(cfg), _g(g) {}

// _____________________________________________________________________________
bool MapConstructor::collapseShrdSegs() { return collapseShrdSegs(DBL_MAX); }

// _____________________________________________________________________________
bool MapConstructor::lineEq(const LineEdge* a, const LineEdge* b) {
  // shortcut
  if (a->pl().getLines().size() != b->pl().getLines().size()) return false;

  const auto shrNd = LineGraph::sharedNode(a, b);

  // TODO: remove quadratic code
  for (const auto& ra : a->pl().getLines()) {
    bool found = false;
    for (const auto& rb : b->pl().getLines()) {
      if (ra.line == rb.line && ra.style == rb.style) {
        if (!shrNd->pl().connOccurs(ra.line, a, b)) return false;

        found = true;

        if (ra.direction == 0 && rb.direction == 0) {
          found = true;
        }
        if (ra.direction == shrNd && rb.direction != 0 &&
            rb.direction != shrNd) {
          found = true;
        }
        if (ra.direction != shrNd && ra.direction != 0 &&
            rb.direction == shrNd) {
          found = true;
        }

        if (found) break;
      }
    }
    if (!found) return false;
  }
  return true;
}

// _____________________________________________________________________________
bool MapConstructor::collapseShrdSegs(double dCut) {
  return collapseShrdSegs(dCut, INT_MAX);
}

// _____________________________________________________________________________
LineNode* MapConstructor::ndCollapseCand(const std::set<LineNode*>& notFrom,
    const size_t numLines,
                                         const double dCut,
                                         const util::geo::Point<double>& point,
                                         const LineNode* spanA,
                                         const LineNode* spanB,
                                         NodeGrid& grid,
                                         LineGraph* g) const {
  LineNode* ndMin = 0;

  std::set<LineNode*> neighbors;

  grid.get(point, dCut * 2, &neighbors);

  double dMin = std::numeric_limits<double>::infinity();

  double dSpanA = std::numeric_limits<double>::infinity();
  double dSpanB = std::numeric_limits<double>::infinity();

  if (spanA) dSpanA = util::geo::dist(point, *spanA->pl().getGeom());
  if (spanB) dSpanB = util::geo::dist(point, *spanB->pl().getGeom());

  for (auto* ndTest : neighbors) {
    if (ndTest->getDeg() == 0) continue;
    if (notFrom.count(ndTest)) continue;
    double d = util::geo::dist(point, *ndTest->pl().getGeom());

    double dMax = maxD(numLines, ndTest, dCut);

    if (d < dSpanA && d < dSpanB && d < dMax && d < dMin) {
      dMin = d;
      ndMin = ndTest;
    }
  }

  LineNode* ret = 0;

  if (ndMin) {
    ndMin->pl().setGeom(util::geo::centroid(
        util::geo::LineSegment<double>(*ndMin->pl().getGeom(), point)));
    grid.remove(ndMin);
    ret = ndMin;
  } else {
    ret = g->addNd(point);
  }

  grid.add(*ret->pl().getGeom(), ret);
  return ret;
}

// _____________________________________________________________________________
double MapConstructor::maxD(size_t lines, const LineNode* nd, double d) const {
    return d;
    size_t numLinesOther = LineGraph::getMaxLineNum(nd);

    return  (d * lines) / 2.0 + (d * numLinesOther) / 2.0;
}

// _____________________________________________________________________________
double MapConstructor::maxD(size_t lines, double d) const {
    return d;
    return  (d * lines);
}

// _____________________________________________________________________________
double MapConstructor::maxD(const LineNode* ndA, const LineNode* ndB, double d) const {
    return d;
    size_t lines = LineGraph::getMaxLineNum(ndA);
    size_t numLinesOther = LineGraph::getMaxLineNum(ndB);

    return (d * lines) / 2.0 + (d * numLinesOther) / 2.0;
}

// _____________________________________________________________________________
bool MapConstructor::collapseShrdSegs(double dCut, size_t steps) {
  // densify

  for (size_t ITER = 0; ITER < 40; ITER++) {
    std::cerr << "ITER " << ITER << std::endl;
    shared::linegraph::LineGraph tgNew;

    // new grid per iteration
    NodeGrid grid(120, 120, bbox());

    std::unordered_map<LineNode*, LineNode*> imageNds;
    std::set<LineNode*> imageNdsSet;

    double SEGL = 10;

    std::vector<std::pair<double, LineEdge*>> sortedEdges;
    for (auto n : *_g->getNds()) {
      for (auto e : n->getAdjList()) {
        if (e->getFrom() != n) continue;
        sortedEdges.push_back({e->pl().getPolyline().getLength(), e});
      }
    }

    std::sort(sortedEdges.rbegin(), sortedEdges.rend());

    size_t j = 0;
    for (const auto& ep : sortedEdges) {
      j++;

      auto e = ep.second;

      LineNode* last = 0;

      std::set<LineNode*> myNds;

      size_t i = 0;
      std::vector<LineNode*> affectedNodes;
      LineNode* front = 0;

      bool imgFromCovered = false;
      bool imgToCovered = false;

      auto pl = *e->pl().getGeom();
      // pl.insert(pl.begin(), *e->getFrom()->pl().getGeom());
      // pl.insert(pl.end(), *e->getTo()->pl().getGeom());

      const auto& plDense = util::geo::densify(util::geo::simplify(*e->pl().getGeom(), 0.5), 10);

      for (const auto& point : plDense) {
        LineNode* cur =  ndCollapseCand(myNds, e->pl().getLines().size(), dCut, point, front, e->getTo(), grid, &tgNew);
        if (i == 0) {
          // this is the "FROM" node
          if (!imageNds.count(e->getFrom())) {
            imageNds[e->getFrom()] = cur;
            imageNdsSet.insert(cur);
            imgFromCovered = true;
          }
        }

        if (i == plDense.size() - 1) {
          // this is the "TO" node
          if (!imageNds.count(e->getTo())) {
            imageNds[e->getTo()] = cur;
            imageNdsSet.insert(cur);
            imgToCovered = true;
          }
        }

        myNds.insert(cur);

        // careful, increase this before the continue below
        i++;

        if (last == cur) continue;  // skip self-edges

        if (cur == imageNds[e->getFrom()]) {
          imgFromCovered = true;
        }
        if (imageNds.count(e->getTo()) && cur == imageNds[e->getTo()]) {
          imgToCovered = true;
        }

        if (last) {
          auto newE = tgNew.getEdg(last, cur);
          if (!newE)  newE = tgNew.addEdg(last, cur);

          combContEdgs(newE, e);
          mergeLines(newE, e, last, cur);
        }

        affectedNodes.push_back(cur);
        if (!front) front = cur;
        last = cur;
      }

      assert(imageNds[e->getFrom()]);
      assert(imageNds[e->getTo()]);

      if (!imgFromCovered) {
        auto newE = tgNew.getEdg(imageNds[e->getFrom()], front);
        if (!newE) newE = tgNew.addEdg(imageNds[e->getFrom()], front);

        combContEdgs(newE, e);
        mergeLines(newE, e, imageNds[e->getFrom()], front);
      }

      if (!imgToCovered) {
        auto newE = tgNew.getEdg(last, imageNds[e->getTo()]);
        if (!newE) newE = tgNew.addEdg(last, imageNds[e->getTo()]);

        combContEdgs(newE, e);
        mergeLines(newE, e, last, imageNds[e->getTo()]);
      }

      // now check all affected nodes for artifact edges (= edges connecting
      // two deg != 1 nodes under the threshold length, they would otherwise
      // never be collapsed because they have to collapse into themself)

      for (const auto& a : affectedNodes) {
        if (imageNdsSet.count(a)) continue;

        double dMin = SEGL;
        LineNode* comb = 0;

        // combine always with the nearest one
        for (auto e : a->getAdjList()) {
          auto b = e->getOtherNd(a);

          if ((a->getDeg() < 3 && b->getDeg() < 3)) continue;
          double dCur = util::geo::dist(*a->pl().getGeom(), *b->pl().getGeom());
          if (dCur <= dMin) {
            dMin = dCur;
            comb = b;
          }
        }

        // this will delete "a" and keep "comb"
        // crucially, "to" has not yet appeared in the list, and we will
        // see the combined node later on
        if (comb && combineNodes(a, comb, &tgNew)) {
          if (a != comb) grid.remove(a);
        }
      }
    }

    std::vector<LineNode*> ndsA;
    ndsA.insert(ndsA.begin(), tgNew.getNds()->begin(), tgNew.getNds()->end());
    for (auto from : ndsA) {
      for (auto e : from->getAdjList()) {
        if (e->getFrom() != from) continue;
        auto to = e->getTo();
        if ((from->getDeg() == 2 || to->getDeg() == 2)) continue;
        if (combineNodes(from, to, &tgNew)) break;
        double dCur =
            util::geo::dist(*from->pl().getGeom(), *to->pl().getGeom());
        if (dCur < maxD(from, to, dCut)) {
          if (combineNodes(from, to, &tgNew)) break;
        }
      }
    }

    // write edge geoms
    for (auto n : *tgNew.getNds()) {
      for (auto e : n->getAdjList()) {
        if (e->getFrom() != n) continue;

        e->pl().setGeom(
            {*e->getFrom()->pl().getGeom(), *e->getTo()->pl().getGeom()});
      }
    }

    // re-collapse
    std::vector<LineNode*> nds;
    nds.insert(nds.begin(), tgNew.getNds()->begin(), tgNew.getNds()->end());

    for (auto n : nds) {
      if (n->getDeg() == 2 &&
          !tgNew.getEdg(n->getAdjList().front()->getOtherNd(n),
                        n->getAdjList().back()->getOtherNd(n))) {
        if (!lineEq(n->getAdjList().front(), n->getAdjList().back())) continue;
        combineEdges(n->getAdjList().front(), n->getAdjList().back(), n,
                     &tgNew);
      }
    }

    // remove edge artifacts
    nds.clear();
    nds.insert(nds.begin(), tgNew.getNds()->begin(), tgNew.getNds()->end());
    for (auto from : nds) {
      for (auto e : from->getAdjList()) {
        if (e->getFrom() != from) continue;
        auto to = e->getTo();
        if (e->pl().getPolyline().getLength() < 2 * maxD(from, to, dCut)) {
          for (auto* oldE : from->getAdjList()) {
            if (e == oldE) continue;

            auto ex = tgNew.getEdg(oldE->getOtherNd(from), to);

            if (ex && ex->pl().getPolyline().getLength() > 2 * maxD(ex->pl().getLines().size(), dCut)) {
                // if long enough, cut the blocking edge in half and add a support node here

                auto plA = ex->pl().getPolyline().getSegment(0, 0.5).getLine();
                auto plB = ex->pl().getPolyline().getSegment(0.5, 1).getLine();
                auto supNd = tgNew.addNd(plA.back());

                auto eA = tgNew.addEdg(ex->getFrom(), supNd, ex->pl());
                auto eB = tgNew.addEdg(supNd, ex->getTo(), ex->pl());

                LineGraph::nodeRpl(eA, ex->getTo(), supNd);
                LineGraph::nodeRpl(eB, ex->getFrom(), supNd);

                eA->pl().setGeom(plA);
                eB->pl().setGeom(plB);

                combContEdgs(eA, ex);
                combContEdgs(eB, ex);

                tgNew.delEdg(ex->getFrom(), ex->getTo());
                delOrigEdgsFor(ex);
            }
          }


          if (combineNodes(from, to, &tgNew)) break;
        }
      }
    }

    // re-collapse again because we might have introduce deg 2 nodes above
    nds.clear();
    nds.insert(nds.begin(), tgNew.getNds()->begin(), tgNew.getNds()->end());

    for (auto n : nds) {
      if (n->getDeg() == 2 &&
          !tgNew.getEdg(n->getAdjList().front()->getOtherNd(n),
                        n->getAdjList().back()->getOtherNd(n))) {
        if (!lineEq(n->getAdjList().front(), n->getAdjList().back())) continue;
        combineEdges(n->getAdjList().front(), n->getAdjList().back(), n,
                     &tgNew);
      }
    }

    // convergence criteria
    double THRESHOLD = 0.001;

    double LEN_OLD = 0;
    double LEN_NEW = 0;
    for (const auto& nd : *_g->getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        LEN_OLD += e->pl().getPolyline().getLength();
      }
    }

    for (const auto& nd : *tgNew.getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        LEN_NEW += e->pl().getPolyline().getLength();
      }
    }

    std::cerr << "len before " << LEN_OLD << std::endl;
    std::cerr << "len after " << LEN_NEW << std::endl;
    std::cerr << "gap " << (1 - LEN_NEW / LEN_OLD) << std::endl;


    if (fabs(1 - LEN_NEW / LEN_OLD) < THRESHOLD) {
      *_g = std::move(tgNew);
      break;
    }

    *_g = std::move(tgNew);
  }

  // std::cerr << "outputting..." << std::endl;

  // // output
  // util::geo::output::GeoGraphJsonOutput out;
  // out.print(*_g, std::cout);

  // exit(0);

  // copy to target graph
  return true;
}

// _____________________________________________________________________________
void MapConstructor::averageNodePositions() {
  for (auto n : *_g->getNds()) {
    double x = 0, y = 0;
    size_t c = 0;

    for (auto e : n->getAdjList()) {
      if (e->getTo() != n) {
        x += e->pl().getPolyline().front().getX();
        y += e->pl().getPolyline().front().getY();
      } else {
        x += e->pl().getPolyline().back().getX();
        y += e->pl().getPolyline().back().getY();
      }
      c++;
    }

    if (c > 0)
      n->pl().setGeom(
          DPoint(x / static_cast<double>(c), y / static_cast<double>(c)));
  }
}

// _____________________________________________________________________________
void MapConstructor::removeEdgeArtifacts() {
  while (contractNodes()) {
  };
}

// _____________________________________________________________________________
void MapConstructor::removeNodeArtifacts() {
  while (contractEdges()) {
  };
}

// _____________________________________________________________________________
bool MapConstructor::contractNodes() {
  for (auto n : *_g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      // contract edges below minimum length, and dead end edges ending in a
      // non-station node
      if (e->pl().getPolyline().getLength() < _cfg->minSegLength) {
        auto from = e->getFrom();
        auto to = e->getTo();
        if (combineNodes(from, to, _g)) return true;
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool MapConstructor::contractEdges() {
  for (auto n : *_g->getNds()) {
    std::vector<LineEdge*> edges;
    edges.insert(edges.end(), n->getAdjList().begin(), n->getAdjList().end());
    if (edges.size() == 2) {
      if (!_g->getEdg(edges[0]->getOtherNd(n), edges[1]->getOtherNd(n))) {
        if (lineEq(edges[0], edges[1])) {
          combineEdges(edges[0], edges[1], n, _g);
          return true;
        }
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool MapConstructor::combineEdges(LineEdge* a, LineEdge* b, LineNode* n) {
  return combineEdges(a, b, n, _g);
}

// _____________________________________________________________________________
bool MapConstructor::combineEdges(LineEdge* a, LineEdge* b, LineNode* n,
                                  LineGraph* g) {
  assert((a->getTo() == n || a->getFrom() == n) &&
         (b->getTo() == n || b->getFrom() == n));

  LineEdge* newEdge = 0;
  util::geo::PolyLine<double> newPl;

  // TODO: there is some copying going on below, which is not always necessary.
  // insert a non-const getLine to polyline and re-use existing polylines where
  // possible

  if (a->getTo() == n && b->getTo() != n) {
    //   a       b
    // ----> n ---->
    auto lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.begin(), lineB.end());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = g->addEdg(a->getFrom(), b->getTo(), a->pl());
    LineGraph::nodeRpl(newEdge, n, newEdge->getTo());
  }

  if (a->getTo() != n && b->getTo() == n) {
    //   a       b
    // <---- n <----
    auto lineB = b->pl().getPolyline().getLine();
    const auto& lineA = a->pl().getPolyline().getLine();
    lineB.insert(lineB.end(), lineA.begin(), lineA.end());
    newPl = util::geo::PolyLine<double>(lineB);

    newEdge = g->addEdg(b->getFrom(), a->getTo(), b->pl());
    LineGraph::nodeRpl(newEdge, n, newEdge->getTo());
  }

  if (a->getFrom() == n && b->getFrom() == n) {
    //   a       b
    // <---- n ---->
    auto lineA = a->pl().getPolyline().getLine();
    std::reverse(lineA.begin(), lineA.end());
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.begin(), lineB.end());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = g->addEdg(a->getTo(), b->getTo(), b->pl());
    LineGraph::nodeRpl(newEdge, n, newEdge->getFrom());
  }

  if (a->getTo() == n && b->getTo() == n) {
    //   a       b
    // ----> n <----
    auto lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.rbegin(), lineB.rend());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = g->addEdg(a->getFrom(), b->getFrom(), a->pl());
    LineGraph::nodeRpl(newEdge, n, newEdge->getTo());
  }

  // set new polyline and simplify a bit
  newPl.simplify(0.5);
  newEdge->pl().setPolyline(newPl);

  combContEdgs(newEdge, a);
  combContEdgs(newEdge, b);

  LineGraph::edgeRpl(a->getFrom(), a, newEdge);
  LineGraph::edgeRpl(a->getTo(), a, newEdge);
  LineGraph::edgeRpl(b->getFrom(), b, newEdge);
  LineGraph::edgeRpl(b->getTo(), b, newEdge);

  delOrigEdgsFor(g->getEdg(a->getFrom(), a->getTo()));
  delOrigEdgsFor(g->getEdg(b->getFrom(), b->getTo()));
  g->delEdg(a->getFrom(), a->getTo());
  g->delEdg(b->getFrom(), b->getTo());

  delOrigEdgsFor(n);
  g->delNd(n);

  return true;
}

// _____________________________________________________________________________
size_t MapConstructor::freeze() {
  size_t i = 0;
  _origEdgs.push_back(OrigEdgs());
  for (auto nd : *_g->getNds()) {
    for (auto* edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      _origEdgs.back()[edg].insert(edg);
      i++;
    }
  }

  return _origEdgs.size() - 1;
}

// _____________________________________________________________________________
void MapConstructor::combContEdgs(const LineEdge* a, const LineEdge* b) {
  for (auto& oe : _origEdgs) {
    oe[a].insert(oe[b].begin(), oe[b].end());
  }
}

// _____________________________________________________________________________
void MapConstructor::delOrigEdgsFor(const LineEdge* a) {
  for (auto& oe : _origEdgs) {
    oe.erase(a);
  }
}

// _____________________________________________________________________________
void MapConstructor::delOrigEdgsFor(const LineNode* a) {
  if (!a) return;
  for (auto* edg : a->getAdjList()) {
    for (auto& oe : _origEdgs) {
      oe.erase(edg);
    }
  }
}

// _____________________________________________________________________________
bool MapConstructor::combineNodes(LineNode* a, LineNode* b) {
  return combineNodes(a, b, _g);
}

// _____________________________________________________________________________
bool MapConstructor::combineNodes(LineNode* a, LineNode* b, LineGraph* g) {
  LineEdge* connecting = g->getEdg(a, b);
  assert(connecting);

  // we will delete a and the connecting edge {a, b}.
  // b will be the contracted node
  //
  b->pl().setGeom(util::geo::centroid(
      util::geo::LineSegment<double>(*a->pl().getGeom(), *b->pl().getGeom())));

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getFrom() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getTo());
    auto* newE = g->getEdg(b, oldE->getTo());

    if (!newE) {
      // add a new edge going from b to the non-a node
      newE = g->addEdg(b, oldE->getTo(), oldE->pl());

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    }

    combContEdgs(newE, oldE);
  }

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getTo() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getFrom());
    auto* newE = g->getEdg(oldE->getFrom(), b);

    if (!newE) {
      newE = g->addEdg(oldE->getFrom(), b, oldE->pl());

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    }

    combContEdgs(newE, oldE);
  }

  delOrigEdgsFor(g->getEdg(a, b));
  g->delEdg(a, b);
  if (a != b) {
    delOrigEdgsFor(a);
    g->delNd(a);
  }

  return true;
}

// _____________________________________________________________________________
PolyLine<double> MapConstructor::geomAvg(const LineEdgePL& geomA, double startA,
                                         double endA, const LineEdgePL& geomB,
                                         double startB, double endB) {
  PolyLine<double> a, b;

  if (startA > endA) {
    a = geomA.getPolyline().getSegment(endA, startA);
    a.reverse();
  } else {
    a = geomA.getPolyline().getSegment(startA, endA);
  }

  if (startB > endB) {
    b = geomB.getPolyline().getSegment(endB, startB);
    b.reverse();
  } else {
    b = geomB.getPolyline().getSegment(startB, endB);
  }

  std::vector<double> weights{
      1.0 * geomA.getLines().size() * geomA.getLines().size(),
      1.0 * geomB.getLines().size() * geomB.getLines().size()};

  PolyLine<double> ret = PolyLine<double>::average({&a, &b});
  ret.simplify(5);
  return ret;
}

// _____________________________________________________________________________
DBox MapConstructor::bbox() const {
  DBox b;

  for (auto n : *_g->getNds()) {
    b = extendBox(*n->pl().getGeom(), b);
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      b = extendBox(e->pl().getPolyline().getLine(), b);
    }
  }

  return b;
}

// _____________________________________________________________________________
bool MapConstructor::foldEdges(LineEdge* a, LineEdge* b) {
  const auto shrNd = LineGraph::sharedNode(a, b);
  assert(shrNd);

  /*
   *
   *                    b
   *           shrNd --------> v
   *            \             /
   *             \           /
   *              \         /
   *             a \       /
   *                \     /
   *                 \   /
   *                  \ /
   *              majNonShrNd
   *
   *
   *   b is the new edge
   */

  if (b->pl().getGeom()->size() == 0 && b->pl().getGeom()->size() == 0) {
    auto v = b->getOtherNd(shrNd);
    v->pl().setGeom(util::geo::centroid(util::geo::LineSegment<double>(
        *v->pl().getGeom(), *a->getOtherNd(shrNd)->pl().getGeom())));
  } else {
    if (b->getTo() == a->getTo() || a->getFrom() == b->getFrom()) {
      b->pl().setPolyline(geomAvg(b->pl(), 0, 1, a->pl(), 0, 1));
    } else {
      b->pl().setPolyline(geomAvg(b->pl(), 0, 1, a->pl(), 1, 0));
    }
  }

  for (auto ro : a->pl().getLines()) {
    if (!b->pl().hasLine(ro.line)) {
      // simply add the route
      if (ro.direction == 0)
        b->pl().addLine(ro.line, 0);
      else if (ro.direction == shrNd)
        b->pl().addLine(ro.line, shrNd);
      else if (ro.direction != shrNd)
        b->pl().addLine(ro.line, b->getOtherNd(shrNd));
    } else {
      auto old = b->pl().lineOcc(ro.line);
      if (ro.direction == 0 && old.direction != 0) {
        // now goes in both directions
        b->pl().delLine(ro.line);
        b->pl().addLine(ro.line, 0);
      }

      if (ro.direction == shrNd && old.direction != shrNd) {
        // now goes in both directions
        b->pl().delLine(ro.line);
        b->pl().addLine(ro.line, 0);
      }

      if (ro.direction != shrNd && old.direction == shrNd) {
        // now goes in both directions
        b->pl().delLine(ro.line);
        b->pl().addLine(ro.line, 0);
      }
    }
  }

  return true;
}

// _____________________________________________________________________________
LineEdgePair MapConstructor::split(LineEdgePL& a, LineNode* fr, LineNode* to,
                                   double p) {
  LineEdge* ret;
  auto right = a.getPolyline().getSegment(p, 1);
  a.setPolyline(a.getPolyline().getSegment(0, p));
  auto helper = _g->addNd(a.getPolyline().back());
  auto helperEdg = _g->addEdg(helper, to, right);

  for (size_t i = 0; i < a.getLines().size(); i++) {
    auto ro = a.getLines()[i];
    if (ro.direction == to) {
      auto* route = ro.line;  // store because of deletion below
      a.delLine(ro.line);
      a.addLine(route, helper);
      helperEdg->pl().addLine(route, to);
      i--;
    } else if (ro.direction == fr) {
      helperEdg->pl().addLine(ro.line, helper);
    } else {
      helperEdg->pl().addLine(ro.line, 0);
    }
  }

  ret = _g->addEdg(fr, helper, a);

  return {ret, helperEdg};
}

// _____________________________________________________________________________
void MapConstructor::mergeLines(LineEdge* newE, LineEdge* oldE, LineNode* newFrom, LineNode* newTo) {
  // add the lines, update the line directions accordingly
  for (auto l : oldE->pl().getLines()) {
    if (!l.direction) {
      newE->pl().addLine(l.line, 0, l.style);
    } else if (l.direction == oldE->getTo()) {
      newE->pl().addLine(l.line, newTo, l.style);
    } else {
      newE->pl().addLine(l.line, newFrom, l.style);
    }
  }
}

// _____________________________________________________________________________
bool MapConstructor::cleanUpGeoms() {
  for (auto n : *_g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      e->pl().setPolyline(e->pl().getPolyline().getSegment(
          e->pl()
              .getPolyline()
              .projectOn(*e->getFrom()->pl().getGeom())
              .totalPos,
          e->pl()
              .getPolyline()
              .projectOn(*e->getTo()->pl().getGeom())
              .totalPos));
    }
  }

  // TODO: edges which continue to each other should be re-connected here

  return true;
}
