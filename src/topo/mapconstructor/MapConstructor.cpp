// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include <climits>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/mapconstructor/MapConstructor.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace topo;

using topo::config::TopoConfig;

using util::geo::Point;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Box;
using util::geo::DBox;
using util::geo::extendBox;
using util::geo::PolyLine;
using util::geo::SharedSegments;

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::Station;
using shared::transitgraph::TransitEdge;
using shared::transitgraph::TransitEdgePair;
using shared::transitgraph::TransitNodePL;
using shared::transitgraph::TransitEdgePL;

double MIN_SEG_LENGTH = 35;
double MAX_SNAP_DIST = 1;

// _____________________________________________________________________________
MapConstructor::MapConstructor(const TopoConfig* cfg, TransitGraph* g)
    : _cfg(cfg), _g(g) {}

// _____________________________________________________________________________
bool MapConstructor::collapseShrdSegs() { return collapseShrdSegs(DBL_MAX); }

// _____________________________________________________________________________
ShrdSegWrap MapConstructor::nextShrdSeg(double dCut, EdgeGrid* grid) {
  double dmin = _cfg->maxAggrDistance;

  for (auto n : *_g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      if (_indEdges.find(e) != _indEdges.end()) continue;

      std::set<TransitEdge*> neighbors;
      grid->getNeighbors(e, fmax(5, dmin * 20), &neighbors);

      for (auto toTest : neighbors) {
        if (_indEdgesPairs.find({e, toTest}) != _indEdgesPairs.end() ||
            _indEdges.find(toTest) != _indEdges.end()) {
          continue;
        }

        if (e != toTest) {
          double dmax =
              fmax(dmin, e->pl().getRoutes().size() * (dmin / 2) +
                             toTest->pl().getRoutes().size() * (dmin / 2));
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
bool MapConstructor::routeEq(const TransitEdge* a, const TransitEdge* b) {
  // shortcut
  if (a->pl().getRoutes().size() != b->pl().getRoutes().size()) return false;

  const auto shrNd = TransitGraph::sharedNode(a, b);

  // TODO: remove quadratic code
  for (const auto& ra : a->pl().getRoutes()) {
    bool found = false;
    for (const auto& rb : b->pl().getRoutes()) {
      if (ra.route == rb.route && ra.style == rb.style) {
        if (!shrNd->pl().connOccurs(ra.route, a, b)) return false;

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
bool MapConstructor::collapseShrdSegs(double dCut, size_t steps) {
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
    TransitNode *a = 0, *b = 0;

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

    TransitEdgePL eaEdgeGeom(ea);
    TransitEdgePL abEdgeGeom(ab);
    TransitEdgePL ecEdgeGeom(ec);
    TransitEdgePL faEdgeGeom(fa);
    TransitEdgePL fcEdgeGeom(fc);

    auto wefrom = w.e->getFrom();
    auto weto = w.e->getTo();
    auto wffrom = w.f->getFrom();
    auto wfto = w.f->getTo();

    for (const auto& r : curEdgeGeom->pl().getRoutes()) {
      if (!r.direction) {
        eaEdgeGeom.addRoute(r.route, 0, r.style);
        abEdgeGeom.addRoute(r.route, 0, r.style);
        ecEdgeGeom.addRoute(r.route, 0, r.style);
      } else if (r.direction == weto) {
        eaEdgeGeom.addRoute(r.route, a, r.style);
        abEdgeGeom.addRoute(r.route, b, r.style);
        ecEdgeGeom.addRoute(r.route, weto, r.style);
      } else {
        eaEdgeGeom.addRoute(r.route, wefrom, r.style);
        abEdgeGeom.addRoute(r.route, a, r.style);
        ecEdgeGeom.addRoute(r.route, b, r.style);
      }
    }

    for (const auto& r : cmpEdgeGeom->pl().getRoutes()) {
      if (!r.direction) {
        faEdgeGeom.addRoute(r.route, 0, r.style);
        abEdgeGeom.addRoute(r.route, 0, r.style);
        fcEdgeGeom.addRoute(r.route, 0, r.style);
      } else if ((r.direction == wfto)) {
        if (fap.totalPos > fbp.totalPos) {
          faEdgeGeom.addRoute(r.route, wfto, r.style);
          abEdgeGeom.addRoute(r.route, a, r.style);
          fcEdgeGeom.addRoute(r.route, b, r.style);
        } else {
          faEdgeGeom.addRoute(r.route, a, r.style);
          abEdgeGeom.addRoute(r.route, b, r.style);
          fcEdgeGeom.addRoute(r.route, wfto, r.style);
        }
      } else {
        if (fap.totalPos > fbp.totalPos) {
          faEdgeGeom.addRoute(r.route, a, r.style);
          abEdgeGeom.addRoute(r.route, b, r.style);
          fcEdgeGeom.addRoute(r.route, wffrom, r.style);
        } else {
          faEdgeGeom.addRoute(r.route, wffrom, r.style);
          abEdgeGeom.addRoute(r.route, a, r.style);
          fcEdgeGeom.addRoute(r.route, b, r.style);
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

    TransitEdge *eaE = 0, *abE = 0, *ebE = 0, *faE = 0, *fbE = 0, *helper = 0;

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
        if ((from->pl().getStops().size() == 0 ||
             to->pl().getStops().size() == 0)) {
          if (combineNodes(from, to)) return true;
        }
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool MapConstructor::contractEdges() {
  for (auto n : *_g->getNds()) {
    std::vector<TransitEdge*> edges;
    edges.insert(edges.end(), n->getAdjList().begin(), n->getAdjList().end());
    if (edges.size() == 2 && n->pl().getStops().size() == 0) {
      if (!_g->getEdg(edges[0]->getOtherNd(n), edges[1]->getOtherNd(n))) {
        if (routeEq(edges[0], edges[1])) {
          combineEdges(edges[0], edges[1], n);
          return true;
        }
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool MapConstructor::combineEdges(TransitEdge* a, TransitEdge* b,
                                  TransitNode* n) {
  assert((a->getTo() == n || a->getFrom() == n) &&
         (b->getTo() == n || b->getFrom() == n));

  TransitEdge* newEdge = 0;
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

    newEdge = _g->addEdg(a->getFrom(), b->getTo(), a->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);
  }

  if (a->getTo() != n && b->getTo() == n) {
    //   a       b
    // <---- n <----
    auto lineB = b->pl().getPolyline().getLine();
    const auto& lineA = a->pl().getPolyline().getLine();
    lineB.insert(lineB.end(), lineA.begin(), lineA.end());
    newPl = util::geo::PolyLine<double>(lineB);

    newEdge = _g->addEdg(b->getFrom(), a->getTo(), b->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);
  }

  if (a->getFrom() == n && b->getFrom() == n) {
    //   a       b
    // <---- n ---->
    auto lineA = a->pl().getPolyline().getLine();
    std::reverse(lineA.begin(), lineA.end());
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.begin(), lineB.end());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = _g->addEdg(a->getTo(), b->getTo(), b->pl());
    routeDirRepl(n, newEdge->getFrom(), newEdge);
  }

  if (a->getTo() == n && b->getTo() == n) {
    //   a       b
    // ----> n <----
    auto lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.rbegin(), lineB.rend());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = _g->addEdg(a->getFrom(), b->getFrom(), a->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);
  }

  // set new polyline and smoothen a bit
  newPl.smoothenOutliers(50);
  newPl.applyChaikinSmooth(3);
  newPl.simplify(3);
  newEdge->pl().setPolyline(newPl);

  combContEdgs(newEdge, a);
  combContEdgs(newEdge, b);

  _g->delEdg(a->getFrom(), a->getTo());
  _g->delEdg(b->getFrom(), b->getTo());
  _g->delNd(n);

  return true;
}

// _____________________________________________________________________________
size_t MapConstructor::freeze() {
  _origEdgs.push_back(OrigEdgs());
  for (auto nd : *_g->getNds()) {
    for (auto* edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      _origEdgs.back()[edg].insert(edg);
    }
  }

  return _origEdgs.size() - 1;
}

// _____________________________________________________________________________
void MapConstructor::combContEdgs(const TransitEdge* a, const TransitEdge* b) {
  for (auto& oe : _origEdgs) {
    oe[a].insert(oe[b].begin(), oe[b].end());
  }
}

// _____________________________________________________________________________
bool MapConstructor::combineNodes(TransitNode* a, TransitNode* b) {
  assert(a->pl().getStops().size() == 0 || b->pl().getStops().size() == 0);

  if (a->pl().getStops().size() != 0) {
    TransitNode* c = a;
    a = b;
    b = c;
  }

  TransitEdge* connecting = _g->getEdg(a, b);
  assert(connecting);

  // we will delete a and the connecting edge {a, b}.
  // b will be the contracted node

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getFrom() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getTo());
    auto* newE = _g->getEdg(b, oldE->getTo());

    if (!newE) {
      // add a new edge going from b to the non-a node
      newE = _g->addEdg(b, oldE->getTo(), oldE->pl());

      // update route dirs
      routeDirRepl(a, b, newE);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      routeDirRepl(a, b, newE);
    }

    combContEdgs(newE, oldE);
  }

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getTo() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getFrom());
    auto* newE = _g->getEdg(oldE->getFrom(), b);

    if (!newE) {
      newE = _g->addEdg(oldE->getFrom(), b, oldE->pl());

      // update route dirs
      routeDirRepl(a, b, newE);
    } else {
      // edge is already existing
      foldEdges(oldE, newE);

      // update route dirs
      routeDirRepl(a, b, newE);
    }

    combContEdgs(newE, oldE);
  }

  _g->delEdg(a, b);
  _g->delNd(a);

  return true;
}

// _____________________________________________________________________________
PolyLine<double> MapConstructor::geomAvg(const TransitEdgePL& geomA,
                                         double startA, double endA,
                                         const TransitEdgePL& geomB,
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
      1.0 * geomA.getRoutes().size() * geomA.getRoutes().size(),
      1.0 * geomB.getRoutes().size() * geomB.getRoutes().size()};

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
void MapConstructor::routeDirRepl(TransitNode* oldN, TransitNode* newN,
                                  TransitEdge* e) {
  auto ro = e->pl().getRoutes().begin();
  while (ro != e->pl().getRoutes().end()) {
    if (ro->direction == oldN) {
      shared::transitgraph::RouteOcc newRo = *ro;
      newRo.direction = newN;

      // delete old
      ro = e->pl().getRoutes().erase(ro);

      // add new
      e->pl().getRoutes().insert(newRo);
    } else {
      ro++;
    }
  }
}

// _____________________________________________________________________________
bool MapConstructor::foldEdges(TransitEdge* a, TransitEdge* b) {
  const auto shrNd = TransitGraph::sharedNode(a, b);
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

  if (b->getTo() == a->getTo() || a->getFrom() == b->getFrom()) {
    b->pl().setPolyline(geomAvg(b->pl(), 0, 1, a->pl(), 0, 1));
  } else {
    b->pl().setPolyline(geomAvg(b->pl(), 0, 1, a->pl(), 1, 0));
  }

  for (auto ro : a->pl().getRoutes()) {
    if (!b->pl().hasRoute(ro.route)) {
      // simply add the route
      if (ro.direction == 0)
        b->pl().addRoute(ro.route, 0);
      else if (ro.direction == shrNd)
        b->pl().addRoute(ro.route, shrNd);
      else if (ro.direction != shrNd)
        b->pl().addRoute(ro.route, b->getOtherNd(shrNd));
    } else {
      auto old = b->pl().getRouteOcc(ro.route);
      if (ro.direction == 0 && old.direction != 0) {
        // now goes in both directions
        b->pl().delRoute(ro.route);
        b->pl().addRoute(ro.route, 0);
      }

      if (ro.direction == shrNd && old.direction != shrNd) {
        // now goes in both directions
        b->pl().delRoute(ro.route);
        b->pl().addRoute(ro.route, 0);
      }

      if (ro.direction != shrNd && old.direction == shrNd) {
        // now goes in both directions
        b->pl().delRoute(ro.route);
        b->pl().addRoute(ro.route, 0);
      }
    }
  }

  return true;
}

// _____________________________________________________________________________
TransitEdgePair MapConstructor::split(TransitEdgePL& a, TransitNode* fr,
                                      TransitNode* to, double p) {
  TransitEdge* ret;
  auto right = a.getPolyline().getSegment(p, 1);
  a.setPolyline(a.getPolyline().getSegment(0, p));
  auto helper = _g->addNd(a.getPolyline().back());
  auto ro = a.getRoutes().begin();
  auto helperEdg = _g->addEdg(helper, to, right);

  while (ro != a.getRoutes().end()) {
    if (ro->direction == to) {
      auto* route = ro->route;  // store because of deletion below
      ro = a.getRoutes().erase(ro);
      a.addRoute(route, helper);
      helperEdg->pl().addRoute(route, to);
    } else if (ro->direction == fr) {
      helperEdg->pl().addRoute(ro->route, helper);
      ro++;
    } else {
      helperEdg->pl().addRoute(ro->route, 0);
      ro++;
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
}

// _____________________________________________________________________________
void MapConstructor::edgeRpl(TransitNode* n, const TransitEdge* oldE,
                             const TransitEdge* newE) {
  if (oldE == newE) return;
  // replace in from
  for (auto& r : n->pl().getConnExc()) {
    auto exFr = r.second.begin();
    while (exFr != r.second.end()) {
      if (exFr->first == oldE) {
        std::swap(r.second[newE], exFr->second);
        exFr = r.second.erase(exFr);
      } else {
        exFr++;
      }
    }
  }

  // replace in to
  for (auto& r : n->pl().getConnExc()) {
    for (auto& exFr : r.second) {
      auto exTo = exFr.second.begin();
      while (exTo != exFr.second.end()) {
        if (*exTo == oldE) {
          exFr.second.insert(newE);
          exTo = exFr.second.erase(exTo);
        } else {
          exTo++;
        }
      }
    }
  }
}
