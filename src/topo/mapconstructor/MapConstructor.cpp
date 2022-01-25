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

double MIN_SEG_LENGTH = 35;
double MAX_SNAP_DIST = 1;

// _____________________________________________________________________________
MapConstructor::MapConstructor(const TopoConfig* cfg, LineGraph* g)
    : _cfg(cfg), _g(g) {}

// _____________________________________________________________________________
bool MapConstructor::collapseShrdSegs() { return collapseShrdSegs(DBL_MAX); }

// _____________________________________________________________________________
ShrdSegWrap MapConstructor::nextShrdSeg(double dCut, EdgeGrid* grid) {
  AggrDistFunc aggrD(_cfg->maxAggrDistance);

  for (auto n : *_g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (_indEdges.find(e) != _indEdges.end()) continue;

      std::set<LineEdge*> neighbors;

      // TODO: why is the * 20 needed here, this does not make sense!
      grid->getNeighbors(e, _cfg->maxAggrDistance * 20, &neighbors);

      for (auto toTest : neighbors) {
        if (_indEdgesPairs.find({e, toTest}) != _indEdgesPairs.end() ||
            _indEdges.find(toTest) != _indEdges.end()) {
          continue;
        }

        if (e != toTest) {
          double dmax = aggrD(e, toTest);
          dmax = fmin(dmax, dCut);

          const auto& s = e->pl().getPolyline().getSharedSegments(
              toTest->pl().getPolyline(), dmax);

          if (s.segments.size() > 0) {
            _pEdges[{e, toTest}]++;
            _pEdges[{toTest, e}]++;

            if (_pEdges[{e, toTest}] > 20) {
              _indEdgesPairs.insert({e, toTest});
              _indEdgesPairs.insert({toTest, e});
            }

            return ShrdSegWrap(e, toTest, s.segments.front());
          }
        }
      }

      // nothing overlapping found for this edge anymore, mark as processed
      _indEdges.insert(e);
    }
  }

  return ShrdSegWrap();
}

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
LineNode* MapConstructor::ndCollapseCand(std::set<LineNode*> notFrom,
                                         double maxD,
                                         const util::geo::Point<double>& point,
                                         const LineNode* spanA,
                                         const LineNode* spanB, NodeGrid& grid,
                                         LineGraph* g) const {
  std::set<LineNode*> neighbors;

  grid.get(point, _cfg->maxAggrDistance * 20, &neighbors);

  double dMin = std::numeric_limits<double>::infinity();
  LineNode* ndMin = 0;
  LineNode* ret = 0;

  double dSpanA = std::numeric_limits<double>::infinity();
  double dSpanB = std::numeric_limits<double>::infinity();

  if (spanA) dSpanA = util::geo::dist(point, *spanA->pl().getGeom());
  if (spanB) dSpanB = util::geo::dist(point, *spanB->pl().getGeom());

  for (auto* ndTest : neighbors) {
    if (ndTest->getDeg() == 0) continue;
    if (notFrom.count(ndTest)) continue;
    double d = util::geo::dist(point, *ndTest->pl().getGeom());
    if (d < dSpanA && d < dSpanB && d < maxD && d < dMin) {
      dMin = d;
      ndMin = ndTest;
    }
  }

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
bool MapConstructor::collapseShrdSegs(double dCut, size_t steps) {
  // densify

  for (size_t ITER = 0; ITER < 40; ITER++) {
    std::cerr << "ITER " << ITER << std::endl;
    shared::linegraph::LineGraph tgNew;
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

    size_t i = 0;
    for (auto ep : sortedEdges) {
      // if (ITER == 11 && i == 36) {
      // // output
      // util::geo::output::GeoGraphJsonOutput out;
      // out.print(*tgNew, std::cout);

      // exit(0);
      // }
      i++;

      auto e = ep.second;
      auto edgePl = e->pl();
      edgePl.setGeom({});

      // TODO change this, making all lines bidirectional here
      for (auto& l : edgePl.getLines()) {
        edgePl.addLine(l.line, 0);
      }

      LineNode* last = 0;

      std::set<LineNode*> myNds;

      auto pl = *e->pl().getGeom();
      pl.insert(pl.begin(), *e->getFrom()->pl().getGeom());
      pl.insert(pl.end(), *e->getTo()->pl().getGeom());

      pl = util::geo::densify(util::geo::simplify(pl, 0.5), 10);
      // pl = util::geo::densify(pl, 10);
      size_t i = 0;
      std::vector<LineNode*> affectedNodes;
      LineNode* front = 0;

      bool imgFromCovered = false;
      bool imgToCovered = false;
      for (const auto& point : pl) {
        LineNode* cur = 0;
        cur =
            ndCollapseCand(myNds, dCut, point, front, e->getTo(), grid, &tgNew);
        if (i == 0) {
          // this is the "FROM" node

          if (!imageNds.count(e->getFrom())) {
            imageNds[e->getFrom()] = cur;
            imageNdsSet.insert(cur);
            imgFromCovered = true;
          }
        }

        if (i == pl.size() - 1) {
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
          auto exE = tgNew.getEdg(last, cur);
          if (exE) {
            for (auto l : edgePl.getLines()) {
              exE->pl().addLine(l.line, l.direction);
            }
            combContEdgs(exE, e);
          } else {
            auto newE = tgNew.addEdg(last, cur, edgePl);
            assert(_origEdgs.back().count(newE) == 0);
            combContEdgs(newE, e);
          }
        }

        affectedNodes.push_back(cur);
        if (!front) front = cur;
        last = cur;
      }

      assert(imageNds[e->getFrom()]);
      assert(imageNds[e->getTo()]);

      if (!imgFromCovered) {
        auto exE = tgNew.getEdg(imageNds[e->getFrom()], front);
        if (exE) {
          for (auto l : edgePl.getLines()) {
            exE->pl().addLine(l.line, l.direction);
          }
          combContEdgs(exE, e);
        } else {
          auto newE = tgNew.addEdg(imageNds[e->getFrom()], front, edgePl);
          combContEdgs(newE, e);
        }
      }

      if (!imgToCovered) {
        auto exE = tgNew.getEdg(last, imageNds[e->getTo()]);
        if (exE) {
          for (auto l : edgePl.getLines()) {
            exE->pl().addLine(l.line, l.direction);
          }
          combContEdgs(exE, e);
        } else {
          auto newE = tgNew.addEdg(last, imageNds[e->getTo()], edgePl);
          combContEdgs(newE, e);
        }
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
        double dCur =
            util::geo::dist(*from->pl().getGeom(), *to->pl().getGeom());
        if (dCur < dCut) {
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

    // output
    // if (ITER == 4) {
    // util::geo::output::GeoGraphJsonOutput out;
    // out.print(*tgNew, std::cout);

    // exit(0);
    // }

    // remove edge artifacts
    nds.clear();
    nds.insert(nds.begin(), tgNew.getNds()->begin(), tgNew.getNds()->end());
    for (auto from : nds) {
      for (auto e : from->getAdjList()) {
        if (e->getFrom() != from) continue;
        if (e->pl().getPolyline().getLength() < dCut) {
          auto to = e->getTo();
          if (combineNodes(from, to, &tgNew)) break;
        }
      }
    }

    std::cerr << "edges before " << _g->numEdgs() << std::endl;
    std::cerr << "edges after " << tgNew.numEdgs() << std::endl;

    size_t THRESHOLD = 2;

    if (fabs(_g->numEdgs() - tgNew.numEdgs()) < THRESHOLD) {
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
bool MapConstructor::collapseShrdSegsOld(double dCut, size_t steps) {
  ShrdSegWrap w;
  _indEdges.clear();
  _indEdgesPairs.clear();
  _pEdges.clear();
  bool found = false;

  EdgeGrid grid = geoIndex();

  while (steps > 0 && (w = nextShrdSeg(dCut, &grid)).e) {
    steps--;
    const auto curEdgeGeom = w.e;
    const auto cmpEdgeGeom = w.f;

    const auto& eap = w.s.first.first;
    const auto& ebp = w.s.second.first;
    const auto& fap = w.s.first.second;
    const auto& fbp = w.s.second.second;

    auto ea = curEdgeGeom->pl().getPolyline().getSegment(0, eap.totalPos);
    auto ec = curEdgeGeom->pl().getPolyline().getSegment(ebp.totalPos, 1);

    PolyLine<double> fa, fc;

    if (fap.totalPos > fbp.totalPos) {
      fa = cmpEdgeGeom->pl().getPolyline().getSegment(fap.totalPos, 1);
      fc = cmpEdgeGeom->pl().getPolyline().getSegment(0, fbp.totalPos);
    } else {
      fa = cmpEdgeGeom->pl().getPolyline().getSegment(0, fap.totalPos);
      fc = cmpEdgeGeom->pl().getPolyline().getSegment(fbp.totalPos, 1);
    }

    auto ab = geomAvg(w.e->pl(), w.s.first.first.totalPos,
                      w.s.second.first.totalPos, w.f->pl(),
                      w.s.first.second.totalPos, w.s.second.second.totalPos);

    // new nodes at the start and end of the shared segment
    LineNode *a = 0, *b = 0;

    if (ea.getLength() < MAX_SNAP_DIST) {
      a = w.e->getFrom();
    }

    if (ec.getLength() < MAX_SNAP_DIST) {
      b = w.e->getTo();
    }

    if (fap.totalPos <= fbp.totalPos) {
      if (fa.getLength() < MAX_SNAP_DIST && b != w.f->getFrom()) {
        a = w.f->getFrom();
      }

      if (fc.getLength() < MAX_SNAP_DIST && a != w.f->getTo()) {
        b = w.f->getTo();
      }
    } else {
      if (fa.getLength() < MAX_SNAP_DIST && b != w.f->getTo()) {
        a = w.f->getTo();
      }

      if (fc.getLength() < MAX_SNAP_DIST && a != w.f->getFrom()) {
        b = w.f->getFrom();
      }
    }

    if (!a) {
      a = _g->addNd({w.s.first.first.p});
    }

    if (!b) {
      b = _g->addNd({w.s.second.first.p});
    }

    if (a == b) continue;

    LineEdgePL eaEdgeGeom(ea);
    LineEdgePL abEdgeGeom(ab);
    LineEdgePL ecEdgeGeom(ec);
    LineEdgePL faEdgeGeom(fa);
    LineEdgePL fcEdgeGeom(fc);

    auto wefrom = w.e->getFrom();
    auto weto = w.e->getTo();
    auto wffrom = w.f->getFrom();
    auto wfto = w.f->getTo();

    for (const auto& r : curEdgeGeom->pl().getLines()) {
      if (!r.direction) {
        eaEdgeGeom.addLine(r.line, 0, r.style);
        abEdgeGeom.addLine(r.line, 0, r.style);
        ecEdgeGeom.addLine(r.line, 0, r.style);
      } else if (r.direction == weto) {
        eaEdgeGeom.addLine(r.line, a, r.style);
        abEdgeGeom.addLine(r.line, b, r.style);
        ecEdgeGeom.addLine(r.line, weto, r.style);
      } else {
        eaEdgeGeom.addLine(r.line, wefrom, r.style);
        abEdgeGeom.addLine(r.line, a, r.style);
        ecEdgeGeom.addLine(r.line, b, r.style);
      }
    }

    for (const auto& r : cmpEdgeGeom->pl().getLines()) {
      if (!r.direction) {
        faEdgeGeom.addLine(r.line, 0, r.style);
        abEdgeGeom.addLine(r.line, 0, r.style);
        fcEdgeGeom.addLine(r.line, 0, r.style);
      } else if ((r.direction == wfto)) {
        if (fap.totalPos > fbp.totalPos) {
          faEdgeGeom.addLine(r.line, wfto, r.style);
          abEdgeGeom.addLine(r.line, a, r.style);
          fcEdgeGeom.addLine(r.line, b, r.style);
        } else {
          faEdgeGeom.addLine(r.line, a, r.style);
          abEdgeGeom.addLine(r.line, b, r.style);
          fcEdgeGeom.addLine(r.line, wfto, r.style);
        }
      } else {
        if (fap.totalPos > fbp.totalPos) {
          faEdgeGeom.addLine(r.line, a, r.style);
          abEdgeGeom.addLine(r.line, b, r.style);
          fcEdgeGeom.addLine(r.line, wffrom, r.style);
        } else {
          faEdgeGeom.addLine(r.line, wffrom, r.style);
          abEdgeGeom.addLine(r.line, a, r.style);
          fcEdgeGeom.addLine(r.line, b, r.style);
        }
      }
    }

    // delete old edges
    grid.remove(_g->getEdg(w.e->getFrom(), w.e->getTo()));
    grid.remove(_g->getEdg(w.f->getFrom(), w.f->getTo()));
    _indEdges.erase(_g->getEdg(w.e->getFrom(), w.e->getTo()));
    _indEdges.erase(_g->getEdg(w.f->getFrom(), w.f->getTo()));
    _g->delEdg(w.e->getFrom(), w.e->getTo());
    _g->delEdg(w.f->getFrom(), w.f->getTo());


    LineEdge *eaE = 0, *abE = 0, *ebE = 0, *faE = 0, *fbE = 0, *helper = 0;

    // add new edges
    if (_g->getEdg(a, b)) {
      std::tie(helper, abE) = split(abEdgeGeom, a, b, .5);
      grid.add(*helper->pl().getGeom(), helper);
      combContEdgs(helper, w.e);
      combContEdgs(helper, w.f);
    } else {
      abE = _g->addEdg(a, b, abEdgeGeom);
    }

    if (a != wefrom) {
      if (_g->getEdg(wefrom, a)) {
        std::tie(helper, eaE) = split(eaEdgeGeom, wefrom, a, .5);
        grid.add(*helper->pl().getGeom(), helper);
        combContEdgs(helper, w.e);
      } else {
        eaE = _g->addEdg(wefrom, a, eaEdgeGeom);
      }
    }

    if (b != weto) {
      if (_g->getEdg(b, weto)) {
        std::tie(ebE, helper) = split(ecEdgeGeom, b, weto, .5);
        grid.add(*helper->pl().getGeom(), helper);
        combContEdgs(helper, w.e);
      } else {
        ebE = _g->addEdg(b, weto, ecEdgeGeom);
      }
    }

    if (fap.totalPos > fbp.totalPos) {
      if (a != wfto) {
        if (_g->getEdg(a, wfto)) {
          std::tie(faE, helper) = split(faEdgeGeom, a, wfto, .5);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          faE = _g->addEdg(a, wfto, faEdgeGeom);
          assert(faE != abE);
        }
      }

      if (b != wffrom) {
        if (_g->getEdg(wffrom, b)) {
          std::tie(helper, fbE) = split(fcEdgeGeom, wffrom, b, .5);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          fbE = _g->addEdg(wffrom, b, fcEdgeGeom);
        }
      }
    } else {
      if (a != wffrom) {
        if (_g->getEdg(wffrom, a)) {
          std::tie(helper, faE) = split(faEdgeGeom, wffrom, a, .5);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          faE = _g->addEdg(wffrom, a, faEdgeGeom);
        }
      }

      if (b != wfto) {
        if (_g->getEdg(b, wfto)) {
          std::tie(fbE, helper) = split(fcEdgeGeom, b, wfto, .5);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          fbE = _g->addEdg(b, wfto, fcEdgeGeom);
        }
      }
    }

    if (abE) {
      grid.add(*abE->pl().getGeom(), abE);
      combContEdgs(abE, w.e);
      combContEdgs(abE, w.f);
    }
    if (eaE) {
      grid.add(*eaE->pl().getGeom(), eaE);
      combContEdgs(eaE, w.e);
    }
    if (ebE) {
      grid.add(*ebE->pl().getGeom(), ebE);
      combContEdgs(ebE, w.e);
    }
    if (faE) {
      grid.add(*faE->pl().getGeom(), faE);
      combContEdgs(faE, w.f);
    }
    if (fbE) {
      grid.add(*fbE->pl().getGeom(), fbE);
      combContEdgs(fbE, w.f);
    }

    found = true;
  }

  return found;
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

  // set new polyline and smoothen a bit
  // newPl.smoothenOutliers(50);
  // newPl.applyChaikinSmooth(1);
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
EdgeGrid MapConstructor::geoIndex() {
  EdgeGrid grid(120, 120, bbox());

  for (auto n : *_g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      grid.add(e->pl().getPolyline().getLine(), e);
    }
  }

  return grid;
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
