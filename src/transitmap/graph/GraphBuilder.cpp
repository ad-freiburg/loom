// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <istream>
#include <set>
#include <stack>
#include <vector>
#include "GraphBuilder.h"
#include "json/json.hpp"
#include "shared/rendergraph/RenderGraph.h"
#include "shared/linegraph/Line.h"
#include "transitmap/config/TransitMapConfig.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using shared::rendergraph::RenderGraph;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::Ordering;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;
using shared::linegraph::NodeFront;
using transitmapper::graph::GraphBuilder;
using util::geo::DPoint;
using util::geo::LinePoint;
using util::geo::LinePointCmp;
using util::geo::PolyLine;

// _____________________________________________________________________________
GraphBuilder::GraphBuilder(const config::Config* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void GraphBuilder::writeMainDirs(RenderGraph* graph) {
  for (auto n : *graph->getNds()) {
    std::set<LineEdge*> eSet;
    eSet.insert(n->getAdjList().begin(), n->getAdjList().end());

    for (LineEdge* e : eSet) {
      NodeFront f(n, e);
      PolyLine<double> pl;

      f.refEtgLengthBefExp = util::geo::len(*e->pl().getGeom());

      if (e->getTo() == n) {
        pl = PolyLine<double>(*e->pl().getGeom())
                 .getOrthoLineAtDist(util::geo::len(*e->pl().getGeom()),
                                     graph->getTotalWidth(e));
      } else {
        pl = PolyLine<double>(*e->pl().getGeom())
                 .getOrthoLineAtDist(0, graph->getTotalWidth(e));
        pl.reverse();
      }

      f.setInitialGeom(pl);

      n->pl().addMainDir(f);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::expandOverlappinFronts(RenderGraph* g) {
  // now, look at the nodes entire front geometries and expand them
  // until nothing overlaps
  double step = 1;

  while (true) {
    bool stillFree = false;
    for (auto n : *g->getNds()) {
      if (n->pl().stops().size() && !g->notCompletelyServed(n) &&
          _cfg->tightStations)
        continue;
      std::set<NodeFront*> overlaps = nodeGetOverlappingFronts(g, n);
      for (auto f : overlaps) {
        stillFree = true;
        if (f->edge->getTo() == n) {
          f->geom = PolyLine<double>(*f->edge->pl().getGeom())
                        .getOrthoLineAtDist(
                            util::geo::len(*f->edge->pl().getGeom()) - step,
                            g->getTotalWidth(f->edge));
        } else {
          f->geom = PolyLine<double>(*f->edge->pl().getGeom())
                        .getOrthoLineAtDist(step, g->getTotalWidth(f->edge));
          f->geom.reverse();
        }

        // cut the edges to fit the new front
        freeNodeFront(n, f);
      }
    }
    if (!stillFree) break;
  }
}

// _____________________________________________________________________________
std::set<NodeFront*> GraphBuilder::nodeGetOverlappingFronts(
    const RenderGraph* g, const LineNode* n) const {
  std::set<NodeFront*> ret;
  double minLength = 6;

  // TODO: why are nodefronts accessed via index?
  for (size_t i = 0; i < n->pl().fronts().size(); ++i) {
    const NodeFront& fa = n->pl().fronts()[i];

    for (size_t j = 0; j < n->pl().fronts().size(); ++j) {
      const NodeFront& fb = n->pl().fronts()[j];

      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      double fac = 1;

      if (nodeFrontsOverlap(g, fa, fb)) {
        if (util::geo::len(*fa.edge->pl().getGeom()) > minLength &&
            fa.geom.distTo(*n->pl().getGeom()) <
                fac * g->getMaxNdFrontWidth(n)) {
          ret.insert(const_cast<NodeFront*>(&fa));
        }
        if (util::geo::len(*fb.edge->pl().getGeom()) > minLength &&
            fb.geom.distTo(*n->pl().getGeom()) <
                fac * g->getMaxNdFrontWidth(n)) {
          ret.insert(const_cast<NodeFront*>(&fb));
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool GraphBuilder::nodeFrontsOverlap(const RenderGraph* g, const NodeFront& a,
                                     const NodeFront& b) const {
  size_t numShr = g->getSharedLines(a.edge, b.edge).size();

  DPoint aa = a.geom.front();
  DPoint ab = a.geom.back();
  DPoint ba = b.geom.front();
  DPoint bb = b.geom.back();

  bool intersects = lineIntersects(aa, ab, ba, bb);
  if (!intersects) return false;

  DPoint i = intersection(aa, ab, ba, bb);

  if (numShr &&
      a.geom.distTo(i) < (g->getWidth(a.edge) + g->getSpacing(a.edge)))
    return true;
  if (b.geom.distTo(a.geom) <
      fmax((g->getWidth(b.edge) + g->getSpacing(b.edge)) * 3,
           (g->getTotalWidth(b.edge) + g->getTotalWidth(a.edge)) / 4))
    return true;

  return false;
}

// _____________________________________________________________________________
void GraphBuilder::freeNodeFront(const LineNode* n, NodeFront* f) {
  PolyLine<double> cutLine = f->geom;

  std::set<LinePoint<double>, LinePointCmp<double>> iSects =
      cutLine.getIntersections(*f->edge->pl().getGeom());
  if (iSects.size() > 0) {
    if (f->edge->getTo() != n) {
      // cut at beginning
      f->edge->pl().setGeom(PolyLine<double>(*f->edge->pl().getGeom())
                                .getSegment(iSects.begin()->totalPos, 1)
                                .getLine());

      assert(cutLine.distTo(f->edge->pl().getGeom()->front()) < 0.1);

    } else {
      // cut at end
      f->edge->pl().setGeom(PolyLine<double>(*f->edge->pl().getGeom())
                                .getSegment(0, (--iSects.end())->totalPos)
                                .getLine());

      assert(cutLine.distTo(f->edge->pl().getGeom()->back()) < 0.1);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::writeInitialConfig(RenderGraph* g) {
  OrderCfg c;
  for (auto n : *g->getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      Ordering order(e->pl().getLines().size());
      for (size_t i = 0; i < e->pl().getLines().size(); i++) order[i] = i;
      c[e] = order;
    }
  }

  g->setConfig(c);
}
