// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include "GraphBuilder.h"
#include "log/Log.h"
#include "./../config/TransitMapConfig.h"

using namespace transitmapper;
using namespace graph;
using namespace gtfsparser;
using namespace gtfs;

using util::geo::Point;

// _____________________________________________________________________________
GraphBuilder::GraphBuilder(const config::Config* cfg)
: _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
}

// _____________________________________________________________________________
void GraphBuilder::writeMainDirs(TransitGraph* graph) {
  for (auto n : *graph->getNodes()) {
    std::set<Edge*> eSet;
    eSet.insert(n->getAdjListIn().begin(), n->getAdjListIn().end());
    eSet.insert(n->getAdjListOut().begin(), n->getAdjListOut().end());

    for (Edge* e : eSet) {
      // atm, always take the first edge trip geometry with the most routes
      size_t lc = 0;

      if (e->getGeom().getLength() == 0) continue;

      NodeFront f(e, n);
      geo::PolyLine pl;

      f.refEtgLengthBefExp = e->getGeom().getLength();

      if (e->getGeomDir() == n) {
        pl = e->getGeom().getOrthoLineAtDist(e->getGeom().getLength(),
            e->getTotalWidth());
      } else {
        pl = e->getGeom().getOrthoLineAtDist(0, e->getTotalWidth());
        pl.reverse();
      }

      f.setGeom(pl);

      // initial free
      freeNodeFront(&f);

      n->addMainDir(f);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::expandOverlappinFronts(TransitGraph* g) {
  // now, look at the nodes entire front geometries and expand them
  // until nothing overlaps
  double step = 1;

  while (true) {
    bool stillFree = false;
    for (auto n : *g->getNodes()) {
      std::set<NodeFront*> overlaps = nodeGetOverlappingFronts(n);
      for (auto f : overlaps) {
        stillFree = true;
        if (f->edge->getGeomDir() == n) {
          f->geom = f->edge->getGeom().getOrthoLineAtDist(
              f->edge->getGeom().getLength() - step,
              f->edge->getTotalWidth());
        } else {
          f->geom = f->edge->getGeom().getOrthoLineAtDist(
              step, f->edge->getTotalWidth());
          f->geom.reverse();
        }

        // cut the edges to fit the new front
        freeNodeFront(f);
      }
    }
    if (!stillFree) break;
  }
}

// _____________________________________________________________________________
std::set<NodeFront*> GraphBuilder::nodeGetOverlappingFronts(const Node* n)
const {
  std::set<NodeFront*> ret;
  double minLength = 1;

  for (size_t i = 0; i < n->getMainDirs().size(); ++i) {
    const NodeFront& fa = n->getMainDirs()[i];

    for (size_t j = 0; j < n->getMainDirs().size(); ++j) {
      const NodeFront& fb = n->getMainDirs()[j];

      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      if ((n->getStops().size() > 0 && fa.geom.distTo(fb.geom) < (fa.edge->getSpacing() + fb.edge->getSpacing()) / 8) ||
          (n->getStops().size() == 0 && nodeFrontsOverlap(fa, fb))) {
        if (fa.edge->getGeom().getLength() > minLength &&
            fa.geom.distTo(n->getPos()) < n->getMaxNodeFrontWidth()) {
          ret.insert(const_cast<NodeFront*>(&fa));
        }
        if (fb.edge->getGeom().getLength() > minLength &&
            fb.geom.distTo(n->getPos()) < n->getMaxNodeFrontWidth()) {
          ret.insert(const_cast<NodeFront*>(&fb));
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool GraphBuilder::nodeFrontsOverlap(const NodeFront& a,
    const NodeFront& b) const {
  size_t numShr= a.edge->getSharedRoutes(*b.edge).size();

  Point aa = a.geom.getLine().front();
  Point ab = a.geom.getLine().back();
  Point ba = b.geom.getLine().front();
  Point bb = b.geom.getLine().back();

  bool intersects = util::geo::lineIntersects(aa, ab, ba, bb);
  if (!intersects) return false;

  Point i = util::geo::intersection(aa, ab, ba, bb);

  if (numShr && a.geom.distTo(i) < (a.edge->getWidth() + a.edge->getSpacing())) return true;
  if (b.geom.distTo(a.geom) < fmax((b.edge->getWidth() + b.edge->getSpacing()) * 3, (b.edge->getTotalWidth() + a.edge->getTotalWidth())/4)) return true;

  return false;
}

// _____________________________________________________________________________
void GraphBuilder::freeNodeFront(NodeFront* f) {
  geo::PolyLine cutLine = f->geom;

  std::set<geo::PointOnLine, geo::PointOnLineCompare> iSects =
      cutLine.getIntersections(f->edge->getGeom());
  if (iSects.size() > 0) {
    if (f->edge->getGeomDir() != f->n) {
      // cut at beginning
      f->edge->setGeom(f->edge->getGeom().getSegment(iSects.begin()->totalPos, 1));
      assert(cutLine.distTo(f->edge->getGeom().getLine().front()) < 0.1);
    } else {
      // cut at end
      f->edge->setGeom(f->edge->getGeom().getSegment(0, (--iSects.end())->totalPos));
      assert(cutLine.distTo(f->edge->getGeom().getLine().back()) < 0.1);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::writeInitialConfig(TransitGraph* g) {
  Configuration c;
  for (graph::Node* n : *g->getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      Ordering order(e->getCardinality());
      for (size_t i = 0; i < e->getCardinality(); i++) order[i] = i;
      c[e] = order;
    }
  }

  g->setConfig(c);
}