// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cassert>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/builder/Builder.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/log/Log.h"

using namespace topo;

using util::geo::Point;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Box;
using util::geo::DBox;
using util::geo::extendBox;
using util::geo::PolyLine;
using util::geo::SharedSegments;

double MIN_SEG_LENGTH = 35;
double MAX_SNAP_DIST = 1;

// _____________________________________________________________________________
Builder::Builder(const config::TopoConfig* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
bool Builder::createTopologicalNodes(TransitGraph* g, bool a) {
  UNUSED(a);
  // TODO: only fallback for tests
  return createTopologicalNodes(g, 999.0);
}

// _____________________________________________________________________________
ShrdSegWrap Builder::nextShrdSeg(TransitGraph* g, double dCut,
                                 EdgeGrid* grid) const {
  double dmin = _cfg->maxAggrDistance;

  for (auto n : *g->getNds()) {
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
bool Builder::routeEq(const TransitEdge* a, const TransitEdge* b) const {
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
bool Builder::createTopologicalNodes(TransitGraph* g, double dCut) {
  return createTopologicalNodes(g, dCut, 10000000);
}

// _____________________________________________________________________________
bool Builder::createTopologicalNodes(TransitGraph* g, double dCut,
                                     size_t steps) {
  ShrdSegWrap w;
  _indEdges.clear();
  _indEdgesPairs.clear();
  _pEdges.clear();
  bool found = false;

  EdgeGrid grid = geoIndex(g);

  while (steps > 0 && (w = nextShrdSeg(g, dCut, &grid)).e) {
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
      a = g->addNd({w.s.first.first.p});
    }

    if (!b) {
      b = g->addNd({w.s.second.first.p});
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
    grid.remove(g->getEdg(w.e->getFrom(), w.e->getTo()));
    grid.remove(g->getEdg(w.f->getFrom(), w.f->getTo()));
    _indEdges.erase(g->getEdg(w.e->getFrom(), w.e->getTo()));
    _indEdges.erase(g->getEdg(w.f->getFrom(), w.f->getTo()));
    g->delEdg(w.e->getFrom(), w.e->getTo());
    g->delEdg(w.f->getFrom(), w.f->getTo());

    TransitEdge *eaE = 0, *abE = 0, *ebE = 0, *faE = 0, *fbE = 0, *helper = 0;

    // add new edges
    if (g->getEdg(a, b)) {
      std::tie(helper, abE) = split(abEdgeGeom, a, b, .5, g);
      grid.add(*helper->pl().getGeom(), helper);
      combContEdgs(helper, w.e);
      combContEdgs(helper, w.f);
    } else {
      abE = g->addEdg(a, b, abEdgeGeom);
    }

    if (a != wefrom) {
      if (g->getEdg(wefrom, a)) {
        std::tie(helper, eaE) = split(eaEdgeGeom, wefrom, a, .5, g);
        grid.add(*helper->pl().getGeom(), helper);
        combContEdgs(helper, w.e);
      } else {
        eaE = g->addEdg(wefrom, a, eaEdgeGeom);
      }
    }

    if (b != weto) {
      if (g->getEdg(b, weto)) {
        std::tie(ebE, helper) = split(ecEdgeGeom, b, weto, .5, g);
        grid.add(*helper->pl().getGeom(), helper);
        combContEdgs(helper, w.e);
      } else {
        ebE = g->addEdg(b, weto, ecEdgeGeom);
      }
    }

    if (fap.totalPos > fbp.totalPos) {
      if (a != wfto) {
        if (g->getEdg(a, wfto)) {
          std::tie(faE, helper) = split(faEdgeGeom, a, wfto, .5, g);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          faE = g->addEdg(a, wfto, faEdgeGeom);
          assert(faE != abE);
        }
      }

      if (b != wffrom) {
        if (g->getEdg(wffrom, b)) {
          std::tie(helper, fbE) = split(fcEdgeGeom, wffrom, b, .5, g);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          fbE = g->addEdg(wffrom, b, fcEdgeGeom);
        }
      }
    } else {
      if (a != wffrom) {
        if (g->getEdg(wffrom, a)) {
          std::tie(helper, faE) = split(faEdgeGeom, wffrom, a, .5, g);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          faE = g->addEdg(wffrom, a, faEdgeGeom);
        }
      }

      if (b != wfto) {
        if (g->getEdg(b, wfto)) {
          std::tie(fbE, helper) = split(fcEdgeGeom, b, wfto, .5, g);
          grid.add(*helper->pl().getGeom(), helper);
          combContEdgs(helper, w.f);
        } else {
          fbE = g->addEdg(b, wfto, fcEdgeGeom);
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
void Builder::averageNodePositions(TransitGraph* g) {
  for (auto n : *g->getNds()) {
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
void Builder::removeEdgeArtifacts(TransitGraph* g) {
  while (contractNodes(g)) {
  };
}

// _____________________________________________________________________________
void Builder::removeNodeArtifacts(TransitGraph* g) {
  while (contractEdges(g)) {
  };
}

// _____________________________________________________________________________
bool Builder::contractNodes(TransitGraph* g) {
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      // contract edges below minimum length, and dead end edges ending in a
      // non-station node
      if (e->pl().getPolyline().getLength() < MIN_SEG_LENGTH) {
        auto from = e->getFrom();
        auto to = e->getTo();
        if ((from->pl().getStops().size() == 0 ||
             to->pl().getStops().size() == 0)) {
          if (combineNodes(from, to, g)) return true;
        }
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool Builder::contractEdges(TransitGraph* g) {
  for (auto n : *g->getNds()) {
    std::vector<TransitEdge*> edges;
    edges.insert(edges.end(), n->getAdjList().begin(), n->getAdjList().end());
    if (edges.size() == 2 && n->pl().getStops().size() == 0) {
      if (!g->getEdg(edges[0]->getOtherNd(n), edges[1]->getOtherNd(n))) {
        if (routeEq(edges[0], edges[1])) {
          combineEdges(edges[0], edges[1], n, g);
          return true;
        }
      }
    }
  }
  return false;
}

// _____________________________________________________________________________
bool Builder::combineEdges(TransitEdge* a, TransitEdge* b, TransitNode* n,
                           TransitGraph* g) {
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

    newEdge = g->addEdg(a->getFrom(), b->getTo(), a->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);
  }

  if (a->getTo() != n && b->getTo() == n) {
    //   a       b
    // <---- n <----
    auto lineB = b->pl().getPolyline().getLine();
    const auto& lineA = a->pl().getPolyline().getLine();
    lineB.insert(lineB.end(), lineA.begin(), lineA.end());
    newPl = util::geo::PolyLine<double>(lineB);

    newEdge = g->addEdg(b->getFrom(), a->getTo(), b->pl());
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

    newEdge = g->addEdg(a->getTo(), b->getTo(), b->pl());
    routeDirRepl(n, newEdge->getFrom(), newEdge);
  }

  if (a->getTo() == n && b->getTo() == n) {
    //   a       b
    // ----> n <----
    auto lineA = a->pl().getPolyline().getLine();
    const auto& lineB = b->pl().getPolyline().getLine();
    lineA.insert(lineA.end(), lineB.rbegin(), lineB.rend());
    newPl = util::geo::PolyLine<double>(lineA);

    newEdge = g->addEdg(a->getFrom(), b->getFrom(), a->pl());
    routeDirRepl(n, newEdge->getTo(), newEdge);
  }

  // set new polyline and smoothen a bit
  newPl.smoothenOutliers(50);
  newPl.applyChaikinSmooth(3);
  newPl.simplify(3);
  newEdge->pl().setPolyline(newPl);

  combContEdgs(newEdge, a);
  combContEdgs(newEdge, b);

  g->delEdg(a->getFrom(), a->getTo());
  g->delEdg(b->getFrom(), b->getTo());
  g->delNd(n);

  return true;
}

// _____________________________________________________________________________
void Builder::initEdges(TransitGraph* g) {
  for (auto nd : *g->getNds()) {
    for (auto* edg : nd->getAdjList()) {
      if (edg->getFrom() != nd) continue;
      _origEdgs[edg].insert(edg);
    }
  }
}

// _____________________________________________________________________________
void Builder::combContEdgs(const TransitEdge* a, const TransitEdge* b) {
  _origEdgs[a].insert(_origEdgs[b].begin(), _origEdgs[b].end());
}

// _____________________________________________________________________________
bool Builder::combineNodes(TransitNode* a, TransitNode* b, TransitGraph* g) {
  assert(a->pl().getStops().size() == 0 || b->pl().getStops().size() == 0);

  if (a->pl().getStops().size() != 0) {
    TransitNode* c = a;
    a = b;
    b = c;
  }

  TransitEdge* connecting = g->getEdg(a, b);
  assert(connecting);

  // we will delete a and the connecting edge {a, b}.
  // b will be the contracted node

  for (auto* oldE : a->getAdjList()) {
    if (oldE->getFrom() != a) continue;
    if (connecting == oldE) continue;

    assert(b != oldE->getTo());
    auto* newE = g->getEdg(b, oldE->getTo());

    if (!newE) {
      // add a new edge going from b to the non-a node
      newE = g->addEdg(b, oldE->getTo(), oldE->pl());

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
    auto* newE = g->getEdg(oldE->getFrom(), b);

    if (!newE) {
      newE = g->addEdg(oldE->getFrom(), b, oldE->pl());

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

  g->delEdg(a, b);
  g->delNd(a);

  return true;
}

// _____________________________________________________________________________
PolyLine<double> Builder::geomAvg(const TransitEdgePL& geomA, double startA,
                                  double endA, const TransitEdgePL& geomB,
                                  double startB, double endB) const {
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
EdgeGrid Builder::geoIndex(const TransitGraph* g) const {
  EdgeGrid grid(120, 120, bbox(g));

  for (auto n : g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      grid.add(e->pl().getPolyline().getLine(), e);
    }
  }

  return grid;
}

// _____________________________________________________________________________
DBox Builder::bbox(const TransitGraph* g) const {
  DBox b;

  for (auto n : g->getNds()) {
    b = extendBox(*n->pl().getGeom(), b);
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      b = extendBox(e->pl().getPolyline().getLine(), b);
    }
  }

  return b;
}

// _____________________________________________________________________________
void Builder::routeDirRepl(TransitNode* oldN, TransitNode* newN,
                           TransitEdge* e) const {
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
bool Builder::foldEdges(TransitEdge* a, TransitEdge* b) {
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
void Builder::collectStations(TransitGraph* g) {
  for (auto nd : *g->getNds()) {
    if (nd->pl().getStops().size()) {
      StationOcc occ{nd->pl().getStops().front(),
                     {nd->getAdjList().begin(), nd->getAdjList().end()}};
      _statClusters.push_back({occ});
      nd->pl().clearStops();
    }
  }
  std::cerr << "Collected " << _statClusters.size() << " station clusters..."
            << std::endl;
}

// _____________________________________________________________________________
std::pair<size_t, size_t> Builder::served(
    const std::vector<TransitEdge*>& adj,
    const std::set<const TransitEdge*>& toServe) {
  std::set<const TransitEdge*> contained;

  for (auto e : adj) contained.insert(_origEdgs[e].begin(), _origEdgs[e].end());

  std::set<const TransitEdge*> iSect;
  set_intersection(contained.begin(), contained.end(), toServe.begin(),
                   toServe.end(), std::inserter(iSect, iSect.begin()));

  // TODO: second return value should be falsely served lines
  return {iSect.size(), 0};
}

// _____________________________________________________________________________
double Builder::candScore(const StationCand& c) {
  double score = 0;

  score += c.dist;
  score += static_cast<double>(c.shouldServ) / static_cast<double>(c.truelyServ) * 100;

  // add a penalty if a station is too close to an existing node
  if (c.edg && c.edg->pl().getPolyline().getLength() * (1 - c.pos) < 100) score += 200;
  if (c.edg && c.edg->pl().getPolyline().getLength() * (c.pos) < 100) score += 200;

  return score;
}

// _____________________________________________________________________________
std::vector<StationCand> Builder::candidates(const StationOcc& occ,
                                             const EdgeGrid& idx) {
  std::vector<StationCand> ret;
  std::set<TransitEdge*> neighbors;
  idx.get(util::geo::pad(util::geo::getBoundingBox(occ.station.pos), 100),
          &neighbors);

  std::cerr << "Got " << neighbors.size() << " candidates..." << std::endl;

  for (auto edg : neighbors) {
    auto pos = edg->pl().getPolyline().projectOn(occ.station.pos);
    double d = util::geo::dist(pos.p, occ.station.pos);

    size_t truelyServed, falselyServed;
    std::tie(truelyServed, falselyServed) = served({edg}, occ.edges);

    ret.push_back(
        StationCand{edg, pos.totalPos, 0, d, occ.edges.size(), truelyServed, falselyServed});

    std::tie(truelyServed, falselyServed) =
        served(edg->getFrom()->getAdjList(), occ.edges);
    ret.push_back(StationCand{
        0, 0, edg->getFrom(),
        util::geo::dist(*edg->getFrom()->pl().getGeom(), occ.station.pos),  occ.edges.size(),
        truelyServed, falselyServed});
    std::tie(truelyServed, falselyServed) =
        served(edg->getTo()->getAdjList(), occ.edges);
    ret.push_back(StationCand{
        0, 0, edg->getTo(),
        util::geo::dist(*edg->getTo()->pl().getGeom(), occ.station.pos),  occ.edges.size(),
        truelyServed, falselyServed});
  }

  struct {
    bool operator()(const StationCand& a, const StationCand& b) const {
      return candScore(a) < candScore(b);
    }
  } cmp;

  std::sort(ret.begin(), ret.end(), cmp);

  std::cerr << "  Cands: " << std::endl;
  for (auto cand : ret) {
    if (cand.edg) {
      std::cerr << "    Edg " << cand.edg << " at position " << cand.pos << " with dist = " << cand.dist << " truely serving " << cand.truelyServ << " edges, falsely serving " << cand.falselyServ << " edges." << std::endl;
    } else {
      std::cerr << "    Nd " << cand.nd << " with dist = " << cand.dist << " truely serving " << cand.truelyServ << " edges, falsely serving " << cand.falselyServ << " edges." << std::endl;
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool Builder::insertStations(TransitGraph* g) {
  auto idx = geoIndex(g);

  for (auto st : _statClusters) {
    if (st.size() == 0) continue;

    auto repr = st.front();
    std::cerr << "Inserting " << repr.station.name << std::endl;

    auto cands = candidates(repr, idx);

    if (cands.size() == 0) continue;

    if (cands.front().edg) {
      auto e = cands.front().edg;

      auto spl = split(e->pl(), e->getFrom(), e->getTo(), cands.front().pos, g);

      combContEdgs(spl.first, e);
      combContEdgs(spl.second, e);

      g->delEdg(e->getFrom(), e->getTo());
      idx.remove(e);
      shared::transitgraph::TransitGraph::sharedNode(spl.first, spl.second)
          ->pl()
          .addStop(repr.station);

      idx.add(*spl.first->pl().getGeom(), spl.first);
      idx.add(*spl.second->pl().getGeom(), spl.second);
    } else {
      cands.front().nd->pl().addStop(repr.station);
    }
  }

  return true;
}

// _____________________________________________________________________________
std::pair<TransitEdge*, TransitEdge*> Builder::split(TransitEdgePL& a,
                                                     TransitNode* fr,
                                                     TransitNode* to, double p,
                                                     TransitGraph* g) const {
  TransitEdge* ret;
  auto right = a.getPolyline().getSegment(p, 1);
  a.setPolyline(a.getPolyline().getSegment(0, p));
  auto helper = g->addNd(a.getPolyline().back());
  auto ro = a.getRoutes().begin();
  auto helperEdg = g->addEdg(helper, to, right);

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

  ret = g->addEdg(fr, helper, a);

  return std::pair<TransitEdge*, TransitEdge*>(ret, helperEdg);
}

// _____________________________________________________________________________
bool Builder::cleanUpGeoms(TransitGraph* g) {
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      e->pl().setPolyline(e->pl().getPolyline().getSegment(e->pl().getPolyline().projectOn(*e->getFrom()->pl().getGeom()).totalPos, e->pl().getPolyline().projectOn(*e->getTo()->pl().getGeom()).totalPos));
    }
  }
}
