// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <istream>
#include <set>
#include <stack>
#include <vector>

#include "GraphBuilder.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;
using shared::linegraph::NodeFront;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::Ordering;
using shared::rendergraph::RenderGraph;
using transitmapper::graph::GraphBuilder;
using util::geo::DPoint;
using util::geo::LinePoint;
using util::geo::LinePointCmp;
using util::geo::PolyLine;

const static double MINL = 10;

// _____________________________________________________________________________
GraphBuilder::GraphBuilder(const config::Config* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void GraphBuilder::writeNodeFronts(RenderGraph* graph) {
  for (auto n : graph->getNds()) {
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

      n->pl().addFront(f);
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::dropOverlappingStations(RenderGraph* graph) {
  util::geo::RTree<LineNode*, util::geo::MultiPolygon, double> tree;

  // store geoms to avoid generating them twice
  std::unordered_map<LineNode*, util::geo::MultiPolygon<double>> geoms;

  // tmp storage of station nodes
  std::vector<LineNode*> stations;

  // fill index
  for (auto n : graph->getNds()) {
    if (!n->pl().stops().size()) continue;
    stations.push_back(n);
    geoms[n] = graph->getStopGeoms(n, _cfg->tightStations, 4);
    tree.add(geoms[n], n);
  }

  double PAD = graph->getWidth(0) + graph->getSpacing(0);

  for (auto n : stations) {
    if (!n->pl().stops().size()) continue;
    std::set<LineNode*> cands;
    auto box = util::geo::pad(util::geo::getBoundingBox(geoms[n]), PAD);
    tree.get(box, &cands);

    for (auto on : cands) {
      if (on == n || on->pl().stops().size() == 0) continue;
      if (util::geo::dist(geoms[n], geoms[on]) <= PAD) {
        // drop n if it is smaller, otherwise wait until the other node
        // is checked
        if (n->getDeg() != 1 &&
            (RenderGraph::getLDeg(n) <= RenderGraph::getLDeg(on) ||
             on->getDeg() == 1)) {
          n->pl().clearStops();
        }
      }
    }
  }
}

// _____________________________________________________________________________
void GraphBuilder::expandOverlappinFronts(RenderGraph* g) {
  // now, look at the nodes entire front geometries and expand them
  // until nothing overlaps
  double step = (g->getWidth(0) + g->getSpacing(0)) / 10;

  while (true) {
    bool stillFree = false;
    for (auto n : g->getNds()) {
      std::set<NodeFront*> overlaps = nodeGetOverlappingFronts(g, n);
      for (auto f : overlaps) {
        stillFree = true;
        double len = util::geo::len(*f->edge->pl().getGeom());

        if (f->edge->getTo() == n) {
          double d = fmax(MINL - 1, len - step);
          f->geom = PolyLine<double>(*f->edge->pl().getGeom())
                        .getOrthoLineAtDist(d, g->getTotalWidth(f->edge));

          f->edge->pl().setGeom(PolyLine<double>(*f->edge->pl().getGeom())
                              .getSegmentAtDist(0, d)
                              .getLine());

        } else {
          double d = fmin(step, len - MINL + 1);
          f->geom = PolyLine<double>(*f->edge->pl().getGeom())
                        .getOrthoLineAtDist(d, g->getTotalWidth(f->edge));

          f->edge->pl().setGeom(PolyLine<double>(*f->edge->pl().getGeom())
                              .getSegmentAtDist(d, len)
                              .getLine());

          f->geom.reverse();
        }

      }
    }
    if (!stillFree) break;
  }
}

// _____________________________________________________________________________
std::set<NodeFront*> GraphBuilder::nodeGetOverlappingFronts(
    const RenderGraph* g, const LineNode* n) const {
  std::set<NodeFront*> ret;

  // TODO: why are nodefronts accessed via index?
  for (size_t i = 0; i < n->pl().fronts().size(); ++i) {
    const NodeFront& fa = n->pl().fronts()[i];

    for (size_t j = i + 1; j < n->pl().fronts().size(); ++j) {
      const NodeFront& fb = n->pl().fronts()[j];

      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      bool overlap = false;

      double maxNfDist = 2 * g->getMaxNdFrontWidth(n);

      if (n->pl().stops().size() && !g->notCompletelyServed(n)) {
        maxNfDist = .5 * g->getMaxNdFrontWidth(n);
        double fac = 0;
        if (_cfg->tightStations) maxNfDist = g->getWidth(0) + g->getSpacing(0);
        overlap = nodeFrontsOverlap(
            g, fa, fb, (g->getWidth(fa.edge) + g->getSpacing(fa.edge)) * fac);
      } else {
        size_t numShr = g->getSharedLines(fa.edge, fb.edge).size();
        double fac = 5;
        if (!numShr) fac = 1;

        overlap = nodeFrontsOverlap(
            g, fa, fb, (g->getWidth(fa.edge) + g->getSpacing(fa.edge)) * fac);
      }

      if (overlap) {
        if (util::geo::len(*fa.edge->pl().getGeom()) > MINL &&
            fa.geom.distTo(*n->pl().getGeom()) < maxNfDist) {
          ret.insert(const_cast<NodeFront*>(&fa));
        }
        if (util::geo::len(*fb.edge->pl().getGeom()) > MINL &&
            fb.geom.distTo(*n->pl().getGeom()) < maxNfDist) {
          ret.insert(const_cast<NodeFront*>(&fb));
        }
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool GraphBuilder::nodeFrontsOverlap(const RenderGraph* g, const NodeFront& a,
                                     const NodeFront& b, double d) const {
  UNUSED(g);
  return b.geom.distTo(a.geom) <= d;
}
