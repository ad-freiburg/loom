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
using shared::linegraph::LineOcc;
using shared::linegraph::Station;

const static double MAX_COLLAPSED_SEG_LENGTH = 500;

// _____________________________________________________________________________
MapConstructor::MapConstructor(const TopoConfig* cfg, LineGraph* g)
    : _cfg(cfg), _g(g) {}

// _____________________________________________________________________________
bool MapConstructor::lineEq(const LineEdge* a, const LineEdge* b) {
  // shortcut
  if (a->pl().getLines().size() != b->pl().getLines().size()) return false;

  const auto shrNd = LineGraph::sharedNode(a, b);

  for (const auto& ra : a->pl().getLines()) {
    if (!b->pl().hasLine(ra.line)) return false;

    const auto& rb = b->pl().lineOcc(ra.line);

    if (!shrNd->pl().connOccurs(ra.line, a, b)) return false;

    bool found = false;

    if (ra.direction == 0 && rb.direction == 0) {
      found = true;
    }
    if (ra.direction == shrNd && rb.direction != 0 && rb.direction != shrNd) {
      found = true;
    }
    if (ra.direction != shrNd && ra.direction != 0 && rb.direction == shrNd) {
      found = true;
    }

    if (!found) return false;
  }
  return true;
}

// _____________________________________________________________________________
LineNode* MapConstructor::ndCollapseCand(
    const std::set<LineNode*>& notFrom, const size_t numLines,
    const double dCut, const util::geo::Point<double>& point,
    const LineNode* spanA, const LineNode* spanB, NodeGeoIdx& geoIdx,
    LineGraph* g) const {
  LineNode* ndMin = 0;

  std::vector<LineNode*> neighbors;

  geoIdx.get(point, dCut, &neighbors);

  double dBest = std::numeric_limits<double>::infinity();

  double dSpanA = std::numeric_limits<double>::infinity();
  double dSpanB = std::numeric_limits<double>::infinity();

  if (spanA)
    dSpanA = util::geo::dist(point, *spanA->pl().getGeom()) / sqrt(2.0);
  if (spanB)
    dSpanB = util::geo::dist(point, *spanB->pl().getGeom()) / sqrt(2.0);

  for (auto* ndTest : neighbors) {
    if (ndTest->getDeg() == 0) continue;
    double d = util::geo::dist(point, *ndTest->pl().getGeom());

    double dMax = maxD(numLines, ndTest, dCut);

    if (d < dSpanA && d < dSpanB && d < dMax && d < dBest &&
        notFrom.count(ndTest) == 0) {
      dBest = d;
      ndMin = ndTest;
    }
  }

  LineNode* ret = 0;

  if (ndMin) {
    ndMin->pl().setGeom(util::geo::centroid(
        util::geo::LineSegment<double>(*ndMin->pl().getGeom(), point)));
    geoIdx.remove(ndMin);
    ret = ndMin;
  } else {
    ret = g->addNd(point);
  }

  geoIdx.add(*ret->pl().getGeom(), ret);
  return ret;
}

// _____________________________________________________________________________
double MapConstructor::maxD(size_t lines, const LineNode* nd, double d) const {
  UNUSED(lines);
  UNUSED(nd);
  return d;
  // size_t numLinesOther = LineGraph::getMaxLineNum(nd);

  // return (d * lines) / 2.0 + (d * numLinesOther) / 2.0;
}

// _____________________________________________________________________________
double MapConstructor::maxD(size_t lines, double d) const {
  UNUSED(lines);
  return d;
  // return (d * lines);
}

// _____________________________________________________________________________
double MapConstructor::maxD(const LineNode* ndA, const LineNode* ndB,
                            double d) const {
  UNUSED(ndA);
  UNUSED(ndB);
  return d;
  // size_t lines = LineGraph::getMaxLineNum(ndA);
  // size_t numLinesOther = LineGraph::getMaxLineNum(ndB);

  // return (d * lines) / 2.0 + (d * numLinesOther) / 2.0;
}

// _____________________________________________________________________________
void MapConstructor::densifyEdg(LineEdge* e, LineGraph* g, double SEGL) {
  return;
  double segd = util::geo::dist(*e->getFrom()->pl().getGeom(),
                                *e->getTo()->pl().getGeom());
  if (segd <= SEGL) return;

  auto a = *e->getFrom()->pl().getGeom();
  auto b = *e->getTo()->pl().getGeom();

  double dx = (a.getX() - b.getX()) / segd;
  double dy = (a.getY() - b.getY()) / segd;
  double curd = SEGL;
  auto last = e->getFrom();
  while (curd < segd) {
    auto p = Point<double>(b.getX() + dx * curd, b.getY() + dy * curd);
    auto nd = g->addNd(p);
    auto newE = g->addEdg(last, nd, e->pl());
    combContEdgs(newE, e);
    LineGraph::nodeRpl(newE, e->getTo(), nd);
    LineGraph::nodeRpl(newE, e->getFrom(), last);
    curd += SEGL;
    last = nd;
  }

  auto newE = g->addEdg(last, e->getTo(), e->pl());
  assert(newE != e);
  combContEdgs(newE, e);
  LineGraph::nodeRpl(newE, e->getFrom(), last);

  g->delEdg(e->getFrom(), e->getTo());
  delOrigEdgsFor(e);
}

// _____________________________________________________________________________
int MapConstructor::collapseShrdSegs() {
  return collapseShrdSegs(_cfg->maxAggrDistance);
}

// _____________________________________________________________________________
int MapConstructor::collapseShrdSegs(double dCut) {
  return collapseShrdSegs(dCut, 50, 5);
}

// _____________________________________________________________________________
int MapConstructor::collapseShrdSegs(double dCut, size_t MAX_ITERS,
                                     double SEGL) {
  size_t ITER = 0;
  for (; ITER < MAX_ITERS; ITER++) {
    shared::linegraph::LineGraph tgNew;

    // new grid per iteration
    NodeGeoIdx geoIdx;

    std::unordered_map<LineNode*, LineNode*> imgNds;
    std::set<LineNode*> imgNdsSet;

    std::vector<std::pair<double, LineEdge*>> sortedEdges;
    for (auto n : _g->getNds()) {
      for (auto e : n->getAdjList()) {
        if (e->getFrom() != n) continue;
        sortedEdges.push_back({e->pl().getPolyline().getLength(), e});
      }
    }

    std::sort(sortedEdges.rbegin(), sortedEdges.rend());

    for (const auto& ep : sortedEdges) {
      auto e = ep.second;

      LineNode* last = 0;

      std::set<LineNode*> myNds;

      size_t i = 0;
      std::vector<LineNode*> affectedNodes;
      LineNode* front = 0;
      LineNode* back = e->getTo();

      bool imgFromCovered = false;
      bool imgToCovered = false;

      util::geo::DLine pl;
      pl.reserve(e->pl().getGeom()->size() + 2);

      pl.push_back(*e->getFrom()->pl().getGeom());
      pl.insert(pl.end(), e->pl().getGeom()->begin(), e->pl().getGeom()->end());
      pl.push_back(*e->getTo()->pl().getGeom());

      const auto& plDense =
          util::geo::densify(util::geo::simplify(pl, 0.5), SEGL);

      for (const auto& point : plDense) {
        if (i == plDense.size() - 1) back = 0;
        LineNode* cur = ndCollapseCand(myNds, e->pl().getLines().size(), dCut,
                                       point, front, back, geoIdx, &tgNew);

        if (i == 0) {
          // this is the "FROM" node
          if (!imgNds.count(e->getFrom())) {
            imgNds[e->getFrom()] = cur;
            imgNdsSet.insert(cur);
            imgFromCovered = true;
          }
        }

        if (i == plDense.size() - 1) {
          // this is the "TO" node
          if (!imgNds.count(e->getTo())) {
            imgNds[e->getTo()] = cur;
            imgNdsSet.insert(cur);
            imgToCovered = true;
          }
        }

        myNds.insert(cur);

        // careful, increase this here, before the continue below
        i++;

        if (last == cur) continue;  // skip self-edges

        if (cur == imgNds[e->getFrom()]) {
          imgFromCovered = true;
        }
        if (imgNds.count(e->getTo()) && cur == imgNds[e->getTo()]) {
          imgToCovered = true;
        }

        if (last) {
          auto newE = tgNew.getEdg(last, cur);
          if (!newE) {
            newE = tgNew.addEdg(last, cur);

            for (const auto& oe : _origEdgs) assert(oe.count(newE) == 0);
          }

          combContEdgs(newE, e);
          mergeLines(newE, e, last, cur);

          densifyEdg(newE, &tgNew, SEGL);
        }

        affectedNodes.push_back(cur);
        if (!front) front = cur;
        last = cur;

        if (imgNds.count(e->getTo()) && last == imgNds.find(e->getTo())->second)
          break;
      }

      assert(imgNds[e->getFrom()]);
      assert(imgNds[e->getTo()]);

      if (!imgFromCovered) {
        auto newE = tgNew.getEdg(imgNds[e->getFrom()], front);
        if (!newE) {
          newE = tgNew.addEdg(imgNds[e->getFrom()], front);

          for (const auto& oe : _origEdgs) assert(oe.count(newE) == 0);
        }

        combContEdgs(newE, e);
        mergeLines(newE, e, imgNds[e->getFrom()], front);

        densifyEdg(newE, &tgNew, SEGL);
      }

      if (!imgToCovered) {
        auto newE = tgNew.getEdg(last, imgNds[e->getTo()]);
        if (!newE) {
          newE = tgNew.addEdg(last, imgNds[e->getTo()]);

          for (const auto& oe : _origEdgs) assert(oe.count(newE) == 0);
        }

        combContEdgs(newE, e);
        mergeLines(newE, e, last, imgNds[e->getTo()]);

        densifyEdg(newE, &tgNew, SEGL);
      }

      // now check all affected nodes for artifact edges (= edges connecting
      // two deg > 2 nodes under the segment length, they would otherwise
      // never be collapsed because they have to collapse into themself)

      for (const auto& a : affectedNodes) {
        if (imgNdsSet.count(a)) continue;

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
        if (comb && combineNodes(a, comb, &tgNew) && a != comb)
          geoIdx.remove(a);
      }
    }

    // soft cleanup
    std::vector<LineNode*> ndsA;
    ndsA.insert(ndsA.begin(), tgNew.getNds().begin(), tgNew.getNds().end());
    for (auto from : ndsA) {
      for (auto e : from->getAdjList()) {
        if (e->getFrom() != from) continue;
        auto to = e->getTo();
        if ((from->getDeg() == 2 || to->getDeg() == 2)) continue;
        if (combineNodes(from, to, &tgNew)) break;
      }
    }

    // write edge geoms
    for (auto n : tgNew.getNds()) {
      for (auto e : n->getAdjList()) {
        if (e->getFrom() != n) continue;

        e->pl().setGeom(
            {*e->getFrom()->pl().getGeom(), *e->getTo()->pl().getGeom()});
      }
    }

    // re-collapse
    std::vector<LineNode*> nds;
    nds.insert(nds.begin(), tgNew.getNds().begin(), tgNew.getNds().end());

    for (auto n : nds) {
      if (n->getDeg() != 2) continue;
      if (!lineEq(n->getAdjList().front(), n->getAdjList().back())) continue;

      // avoid edges that are too long, this will cause
      // orig edge creep!
      if (util::geo::longerThan(*n->getAdjList().front()->pl().getGeom(),
                                *n->getAdjList().back()->pl().getGeom(),
                                MAX_COLLAPSED_SEG_LENGTH))
        continue;

      auto ex = tgNew.getEdg(n->getAdjList().front()->getOtherNd(n),
                             n->getAdjList().back()->getOtherNd(n));

      if (ex && ex->pl().getPolyline().longerThan(
                    2 * maxD(ex->pl().getLines().size(), dCut))) {
        // if long enough, cut the blocking edge in half and add a support
        // node here

        supportEdge(ex, &tgNew);
      } else if (ex) {
        // else dont contract
        continue;
      }

      combineEdges(n->getAdjList().front(), n->getAdjList().back(), n, &tgNew);
    }

    // remove edge artifacts as long as possible, because an artifact removal
    // might introduce another artifact if we fold edges
    bool found;
    do {
      found = false;
      nds.clear();
      nds.insert(nds.begin(), tgNew.getNds().begin(), tgNew.getNds().end());
      std::set<const LineNode*> skip;

      for (auto from : nds) {
        if (skip.count(from)) continue;
        for (auto e : from->getAdjList()) {
          if (e->getFrom() != from) continue;

          auto to = e->getTo();

          if (util::geo::dist(*from->pl().getGeom(), *to->pl().getGeom()) <
                  maxD(from, to, dCut) &&
              e->pl().getPolyline().shorterThan(maxD(from, to, dCut))) {
            for (auto* oldE : from->getAdjList()) {
              if (e == oldE) continue;
              auto ex = tgNew.getEdg(oldE->getOtherNd(from), to);

              if (ex && ex->pl().getPolyline().longerThan(
                            2 * maxD(ex->pl().getLines().size(), dCut))) {
                // if long enough, cut the blocking edge in half and add a
                // support node here
                supportEdge(ex, &tgNew);
              }
            }

            // don't contract cases where
            //   1) one of the nodes is a terminus, the other is degree 2,
            //      and the lines on both edges match
            //   2) both nodes are degree 2 nodes, and the lines on the three
            //      edges match
            // this is to avoid deleting edges that are left during the node 2
            // contraction phase above to avoid too long edges (resulting in
            // line creep)
            if (from->getDeg() == 1 && to->getDeg() == 2) {
              if (lineEq(to->getAdjList().front(), to->getAdjList().back()))
                continue;
            } else if (from->getDeg() == 2 && to->getDeg() == 1) {
              if (lineEq(from->getAdjList().front(), from->getAdjList().back()))
                continue;
              continue;
            } else if (from->getDeg() == 2 && to->getDeg() == 2) {
              if (lineEq(from->getAdjList().front(),
                         from->getAdjList().back()) &&
                  lineEq(to->getAdjList().front(), to->getAdjList().back()))
                continue;
            }

            if (combineNodes(from, to, &tgNew)) {
              found = true;
              break;
            }
          }
        }
      }
    } while (found);

    // re-collapse again because we might have introduced deg 2 nodes above
    nds.clear();
    nds.insert(nds.begin(), tgNew.getNds().begin(), tgNew.getNds().end());

    for (auto n : nds) {
      if (n->getDeg() == 2 &&
          !tgNew.getEdg(n->getAdjList().front()->getOtherNd(n),
                        n->getAdjList().back()->getOtherNd(n))) {
        // avoid edges that are too long, this will cause
        // orig edge creep
        if (util::geo::longerThan(*n->getAdjList().front()->pl().getGeom(),
                                  *n->getAdjList().back()->pl().getGeom(),
                                  MAX_COLLAPSED_SEG_LENGTH))
          continue;

        if (!lineEq(n->getAdjList().front(), n->getAdjList().back())) continue;
        combineEdges(n->getAdjList().front(), n->getAdjList().back(), n,
                     &tgNew);
      }
    }

    // smoothen a bit
    for (auto n : tgNew.getNds()) {
      for (auto e : n->getAdjList()) {
        if (e->getFrom() != n) continue;
        auto& pl = e->pl().getPolyline();
        pl.smoothenOutliers(50);
        pl.simplify(1);
        pl = PolyLine<double>(util::geo::densify(pl.getLine(), 5));
        pl.applyChaikinSmooth(1);
        pl.simplify(1);
      }
    }

    // convergence criteria
    double THRESHOLD = 0.002;

    double LEN_OLD = 0;
    double LEN_NEW = 0;
    for (const auto& nd : _g->getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        LEN_OLD += e->pl().getPolyline().getLength();
      }
    }

    for (const auto& nd : tgNew.getNds()) {
      for (const auto& e : nd->getAdjList()) {
        if (e->getFrom() != nd) continue;
        LEN_NEW += e->pl().getPolyline().getLength();
      }
    }

    *_g = std::move(tgNew);

    LOGTO(DEBUG, std::cerr)
        << "iter " << ITER << ", distance gap: " << (1 - LEN_NEW / LEN_OLD);
    if (fabs(1 - LEN_NEW / LEN_OLD) < THRESHOLD) break;
  }

  return ITER + 1;
}

// _____________________________________________________________________________
void MapConstructor::averageNodePositions() {
  for (auto n : _g->getNds()) {
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
void MapConstructor::removeNodeArtifacts(bool keepStations) {
  while (contractEdges(keepStations)) {
  };
}

// _____________________________________________________________________________
bool MapConstructor::contractNodes() {
  for (auto n : _g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      // contract edges below minimum length
      if (e->pl().getPolyline().shorterThan(_cfg->maxAggrDistance)) {
        auto from = e->getFrom();
        auto to = e->getTo();

        bool dontContract = false;

        // check if we would fold edges with vastly different geoms
        for (auto* oldE : from->getAdjList()) {
          if (e == oldE) continue;

          auto* newE = _g->getEdg(to, oldE->getTo());

          if (newE && fabs(util::geo::len(*newE->pl().getGeom()) -
                           util::geo::len(*oldE->pl().getGeom())) >
                          _cfg->maxAggrDistance * 2) {
            dontContract = true;
          }
        }

        if (!dontContract && combineNodes(from, to, _g)) return true;
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool MapConstructor::contractEdges(bool keepStations) {
  for (auto n : _g->getNds()) {
    if (keepStations && n->pl().stops().size()) continue;
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

  if (a->getTo() == n && b->getTo() != n) {
    //   a       b
    // ----> n ---->
    auto& lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.begin(), lineB.end());
    a->pl().getPolyline().simplify(0.5);

    newEdge = g->addEdg(a->getFrom(), b->getTo(), std::move(a->pl()));
    LineGraph::nodeRpl(newEdge, n, newEdge->getTo());
  }

  if (a->getTo() != n && b->getTo() == n) {
    //   a       b
    // <---- n <----
    auto& lineB = b->pl().getPolyline().getLine();
    const auto& lineA = a->pl().getPolyline().getLine();
    lineB.insert(lineB.end(), lineA.begin(), lineA.end());
    b->pl().getPolyline().simplify(0.5);

    newEdge = g->addEdg(b->getFrom(), a->getTo(), std::move(b->pl()));
    LineGraph::nodeRpl(newEdge, n, newEdge->getTo());
  }

  if (a->getFrom() == n && b->getFrom() == n) {
    //   a       b
    // <---- n ---->
    auto& lineB = b->pl().getPolyline().getLine();
    const auto& lineA = a->pl().getPolyline().getLine();
    std::reverse(lineB.begin(), lineB.end());

    lineB.insert(lineB.end(), lineA.begin(), lineA.end());
    std::reverse(lineB.begin(), lineB.end());
    b->pl().getPolyline().simplify(0.5);

    newEdge = g->addEdg(a->getTo(), b->getTo(), std::move(b->pl()));
    LineGraph::nodeRpl(newEdge, n, newEdge->getFrom());
  }

  if (a->getTo() == n && b->getTo() == n) {
    //   a       b
    // ----> n <----
    auto& lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.rbegin(), lineB.rend());
    a->pl().getPolyline().simplify(0.5);

    newEdge = g->addEdg(a->getFrom(), b->getFrom(), std::move(a->pl()));
    LineGraph::nodeRpl(newEdge, n, newEdge->getTo());
  }

  assert(newEdge != a);
  assert(newEdge != b);

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

  assert(newEdge->getFrom() != n && newEdge->getTo() != n);

  delOrigEdgsFor(n);
  g->delNd(n);

  return true;
}

// _____________________________________________________________________________
size_t MapConstructor::freeze() {
  _origEdgs.push_back(OrigEdgs());
  for (auto nd : _g->getNds()) {
    for (auto* edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      _origEdgs.back()[edg].insert(edg);
    }
  }

  return _origEdgs.size() - 1;
}

// _____________________________________________________________________________
void MapConstructor::combContEdgs(const LineEdge* a, const LineEdge* b) {
  for (auto& oe : _origEdgs) oe[a].insert(oe[b].begin(), oe[b].end());
}

// _____________________________________________________________________________
void MapConstructor::delOrigEdgsFor(const LineEdge* a) {
  for (auto& oe : _origEdgs) oe.erase(a);
}

// _____________________________________________________________________________
void MapConstructor::delOrigEdgsFor(const LineNode* a) {
  if (!a) return;
  for (auto* edg : a->getAdjList()) delOrigEdgsFor(edg);
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

  b->pl().setGeom(util::geo::centroid(
      util::geo::LineSegment<double>(*a->pl().getGeom(), *b->pl().getGeom())));

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getFrom() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getTo());
    auto* newE = g->getEdg(b, oldE->getTo());

    if (!newE) {
      // add a new edge going from b to the non-b node
      newE = g->addEdg(b, oldE->getTo(), std::move(oldE->pl()));

      for (const auto& oe : _origEdgs) assert(oe.count(newE) == 0);

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    }

    assert(newE != connecting);

    combContEdgs(newE, oldE);
    combContEdgs(newE, connecting);
  }

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getTo() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getFrom());
    auto* newE = g->getEdg(oldE->getFrom(), b);

    if (!newE) {
      newE = g->addEdg(oldE->getFrom(), b, std::move(oldE->pl()));

      for (const auto& oe : _origEdgs) assert(oe.count(newE) == 0);

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      LineGraph::nodeRpl(newE, a, b);
    }

    assert(newE != connecting);

    combContEdgs(newE, oldE);
    combContEdgs(newE, connecting);
  }

  // also add to the edges on b, for symmetry reasons
  for (auto* oldE : b->getAdjList()) {
    combContEdgs(oldE, connecting);
  }

  assert(g->getEdg(a, b));
  assert(g->getEdg(a, b) == connecting);

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
  ret.simplify(0.5);
  return ret;
}

// _____________________________________________________________________________
DBox MapConstructor::bbox() const {
  DBox b;

  for (auto n : _g->getNds()) {
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
      } else if (ro.direction == shrNd && old.direction != shrNd) {
        // now goes in both directions
        b->pl().delLine(ro.line);
        b->pl().addLine(ro.line, 0);
      } else if (ro.direction != shrNd && old.direction == shrNd) {
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
void MapConstructor::mergeLines(LineEdge* newE, LineEdge* oldE,
                                LineNode* newFrom, LineNode* newTo) {
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
  for (auto n : _g->getNds()) {
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

// _____________________________________________________________________________
void MapConstructor::removeOrphanLines() {
  bool found;

  do {
    found = false;
    std::vector<LineEdge*> toDelEdgs;

    for (auto n : _g->getNds()) {
      for (auto e : n->getAdjList()) {
        if (e->getFrom() != n) continue;

        std::vector<const shared::linegraph::Line*> toDel;

        for (const auto& lo : e->pl().getLines()) {
          if (((e->getFrom()->pl().stops().size() == 0 ||
                !e->getFrom()->pl().lineServed(lo.line)) &&
               LineGraph::terminatesAt(e, e->getFrom(), lo.line)) ||
              ((e->getTo()->pl().stops().size() == 0 ||
                !e->getTo()->pl().lineServed(lo.line)) &&
               LineGraph::terminatesAt(e, e->getTo(), lo.line))) {
            toDel.push_back(lo.line);
          }
        }

        // if the edge would run empty, and isnt a stump, dont delete
        if (e->getFrom()->getDeg() != 1 && e->getTo()->getDeg() != 1 &&
            toDel.size() == e->pl().getLines().size()) {
          continue;
        }

        // check if this edge was used as a path to connect two edges with this
        // line
        for (const auto& del : toDel) {
          for (auto e2 : n->getAdjList()) {
            if (e == e2) continue;
            for (auto e3 : n->getAdjList()) {
              if (e3 == e2) continue;
              if (LineGraph::lineCtd(e, e2, del) &&
                  LineGraph::lineCtd(e, e3, del) &&
                  !LineGraph::lineCtd(e2, e3, del)) {
                n->pl().delConnExc(del, e2, e3);
              }
            }
          }

          // clear restrictions
          for (auto other : e->getFrom()->getAdjList())
            e->getFrom()->pl().delConnExc(del, e, other);
          for (auto other : e->getTo()->getAdjList())
            e->getTo()->pl().delConnExc(del, e, other);

          e->pl().delLine(del);
          found = true;
        }

        // if the edge runs empty, delete it
        if (e->pl().getLines().size() == 0) toDelEdgs.push_back(e);
      }

      // del removed lines which were marked served
      std::vector<const shared::linegraph::Line*> notServedToDel;
      std::set<const shared::linegraph::Line*> adjLines;
      for (auto e : n->getAdjList()) {
        for (auto l : e->pl().getLines()) {
          adjLines.insert(l.line);
        }
      }
      for (auto notServed : n->pl().getLinesNotServed()) {
        if (!adjLines.count(notServed)) notServedToDel.push_back(notServed);
      }

      for (auto l : notServedToDel) n->pl().delLineNotServed(l);
    }

    for (auto e : toDelEdgs) _g->delEdg(e->getFrom(), e->getTo());

    std::vector<LineNode*> toDelNds;
    for (auto nd : _g->getNds()) {
      if (nd->getDeg() == 0) toDelNds.push_back(nd);
    }

    for (auto nd : toDelNds) _g->delNd(nd);
  } while (found);
}

// _____________________________________________________________________________
void MapConstructor::reconstructIntersections() {
  averageNodePositions();
  for (auto n : _g->getNds()) {
    for (auto e : n->getAdjList()) {
      auto pl = e->pl().getPolyline();
      pl = pl.getSegmentAtDist(_cfg->maxAggrDistance,
                               pl.getLength() - _cfg->maxAggrDistance);
      auto l = pl.getLine();
      l.insert(l.begin(), *e->getFrom()->pl().getGeom());
      l.insert(l.end(), *e->getTo()->pl().getGeom());
      e->pl().setGeom(l);
    }
  }
}

// _____________________________________________________________________________
void MapConstructor::supportEdge(LineEdge* ex, LineGraph* g) {
  // if long enough, cut the blocking edge in half and add a support
  // node here

  auto plA = ex->pl().getPolyline().getSegment(0, 0.5).getLine();
  auto plB = ex->pl().getPolyline().getSegment(0.5, 1).getLine();
  auto supNd = g->addNd(plA.back());

  auto eA = g->addEdg(ex->getFrom(), supNd, ex->pl());
  auto eB = g->addEdg(supNd, ex->getTo(), ex->pl());

  for (const auto& oe : _origEdgs) assert(oe.count(eB) == 0);
  for (const auto& oe : _origEdgs) assert(oe.count(eA) == 0);

  assert(eA != ex);
  assert(eB != ex);

  combContEdgs(eA, ex);
  combContEdgs(eB, ex);

  LineGraph::nodeRpl(eA, ex->getTo(), supNd);
  LineGraph::nodeRpl(eB, ex->getFrom(), supNd);

  eA->pl().setGeom(plA);
  eB->pl().setGeom(plB);

  g->delEdg(ex->getFrom(), ex->getTo());
  delOrigEdgsFor(ex);
}
