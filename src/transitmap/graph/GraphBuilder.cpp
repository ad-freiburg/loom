// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <proj_api.h>
#include <istream>
#include <set>
#include <stack>
#include <vector>
#include "GraphBuilder.h"
#include "json/json.hpp"
#include "shared/linegraph/Route.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/RenderGraph.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace util::geo;
using transitmapper::graph::GraphBuilder;
using shared::linegraph::NodeFront;
using shared::linegraph::LineNode;
using shared::linegraph::LineEdge;
using shared::linegraph::Route;
using util::geo::PolyLine;

const static char* WGS84_PROJ =
    "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs";

// _____________________________________________________________________________
GraphBuilder::GraphBuilder(const config::Config* cfg) : _cfg(cfg) {
  _mercProj = pj_init_plus(WGS84_PROJ);
}

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
      if (n->pl().getStops().size() && _cfg->tightStations) continue;
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
void GraphBuilder::createMetaNodes(RenderGraph* g) {
  std::vector<NodeFront> cands;
  while ((cands = getNextMetaNodeCand(g)).size() > 0) {
    // remove all edges completely contained
    for (auto nf : cands) {
      auto onfs = getClosedNodeFronts(g, nf.n);

      for (auto onf : onfs) {
        const LineEdge* e = onf.edge;

        bool found = false;

        for (auto nf : cands) {
          if (nf.n == e->getFrom()) {
            found = true;
            break;
          }
        }

        if (!found) continue;

        for (auto nf : cands) {
          if (nf.n == e->getTo()) {
            found = true;
            break;
          }
        }

        if (!found) continue;

        e->getTo()->pl().delMainDir(e);
        e->getFrom()->pl().delMainDir(e);
        g->delEdg(e->getTo(), e->getFrom());
      }
    }

    // first node has new ref node id
    LineNode* ref = g->addNd(*cands[0].n->pl().getGeom());

    std::set<LineNode*> toDel;

    for (auto nf : cands) {
      for (auto onf : getOpenNodeFronts(g, nf.n)) {
        LineEdge* e;
        LineNode* other;
        if (onf.edge->getTo() == nf.n) {
          e = g->addEdg(onf.edge->getFrom(), ref, onf.edge->pl());
          other = onf.edge->getFrom();
        } else {
          e = g->addEdg(ref, onf.edge->getTo(), onf.edge->pl());
          other = onf.edge->getTo();
        }

        // update the directions, if necessary
        std::set<shared::linegraph::RouteOcc> del;
        for (auto& to : e->pl().getRoutes()) {
          if (to.direction == nf.n) del.insert(to);
        }

        // also update the edge of the other node front
        NodeFront* otherFr = other->pl().getNodeFrontFor(onf.edge);
        assert(otherFr);
        otherFr->edge = e;

        // remove the original edge
        g->delEdg(onf.edge->getFrom(), onf.edge->getTo());

        // update the original edge in the checked node front
        onf.edge = e;

        // update the routes
        for (auto ro : del) {
          e->pl().delRoute(ro.route);
          e->pl().addRoute(ro.route, ref);
        }

        toDel.insert(onf.n);

        // update the node front node
        onf.n = ref;

        // add the new node front to the new node
        ref->pl().addMainDir(onf);
      }
    }

    // delete the nodes marked for deletion
    for (auto delNd : toDel) g->delNd(delNd);
  }
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getNextMetaNodeCand(RenderGraph* g) const {
  for (auto n : *g->getNds()) {
    if (n->pl().getStops().size()) continue;
    if (getOpenNodeFronts(g, n).size() != 1) continue;

    std::set<const LineNode*> potClique;

    std::stack<const LineNode*> nodeStack;
    nodeStack.push(n);

    while (!nodeStack.empty()) {
      const LineNode* n = nodeStack.top();
      nodeStack.pop();

      if (n->pl().getStops().size() == 0) {
        potClique.insert(n);
        for (auto nff : getClosedNodeFronts(g, n)) {
          const LineNode* m;

          if (nff.edge->getTo() == n) {
            m = nff.edge->getFrom();
          } else {
            m = nff.edge->getTo();
          }

          if (potClique.find(m) == potClique.end()) {
            nodeStack.push(m);
          }
        }
      }
    }

    if (isClique(g, potClique)) {
      std::vector<NodeFront> ret;

      for (auto n : potClique) {
        if (getOpenNodeFronts(g, n).size() > 0) {
          ret.push_back(getOpenNodeFronts(g, n)[0]);
        } else {
          for (auto nf : getClosedNodeFronts(g, n)) {
            ret.push_back(nf);
          }
        }
      }

      return ret;
    }
  }

  return std::vector<NodeFront>();
}

// _____________________________________________________________________________
bool GraphBuilder::isClique(const RenderGraph* g,
                            std::set<const LineNode*> potClique) const {
  if (potClique.size() < 2) return false;

  for (const LineNode* a : potClique) {
    for (const LineNode* b : potClique) {
      if (util::geo::dist(*a->pl().getGeom(), *b->pl().getGeom()) >
          (_cfg->lineWidth + _cfg->lineSpacing) * 10) {
        return false;
      }
    }
  }

  std::set<const LineNode*> periphery;

  for (const LineNode* n : potClique) {
    for (auto nf : getClosedNodeFronts(g, n)) {
      if (nf.edge->getTo() == n) {
        if (potClique.find(nf.edge->getFrom()) == potClique.end()) {
          return false;
        }
      } else {
        if (potClique.find(nf.edge->getTo()) == potClique.end()) {
          return false;
        }
      }
    }

    // catch cases where two clique nodes share the same node at teir
    // open front
    for (auto nf : getOpenNodeFronts(g, n)) {
      if (nf.edge->getTo() == n) {
        if (periphery.count(nf.edge->getFrom())) return false;
        periphery.insert(nf.edge->getFrom());
      } else {
        if (periphery.count(nf.edge->getTo())) return false;
        periphery.insert(nf.edge->getTo());
      }
    }

    for (auto nf : getOpenNodeFronts(g, n)) {
      if (nf.edge->getTo() == n) {
        if (potClique.find(nf.edge->getFrom()) != potClique.end()) {
          return false;
        }
      } else {
        if (potClique.find(nf.edge->getTo()) != potClique.end()) {
          return false;
        }
      }
    }
  }

  return true;
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getOpenNodeFronts(
    const RenderGraph* g, const LineNode* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->pl().getMainDirs()) {
    if (util::geo::len(*nf.edge->pl().getGeom()) >
            (g->getWidth(nf.edge) + g->getSpacing(nf.edge)) ||
        (nf.edge->getOtherNd(n)->pl().getNodeFrontFor(nf.edge)->geom.distTo(
             *nf.edge->getOtherNd(n)->pl().getGeom()) >
         6 * (g->getWidth(nf.edge) + g->getSpacing(nf.edge))) ||
        (nf.edge->getTo()->pl().getStops().size() > 0) ||
        (nf.edge->getFrom()->pl().getStops().size() > 0)) {
      res.push_back(nf);
    }
  }

  return res;
}

// _____________________________________________________________________________
std::vector<NodeFront> GraphBuilder::getClosedNodeFronts(
    const RenderGraph* g, const LineNode* n) const {
  std::vector<NodeFront> res;
  for (auto nf : n->pl().getMainDirs()) {
    if (!(util::geo::len(*nf.edge->pl().getGeom()) >
          (g->getWidth(nf.edge) + g->getSpacing(nf.edge))) &&
        !(nf.edge->getOtherNd(n)->pl().getNodeFrontFor(nf.edge)->geom.distTo(
              *nf.edge->getOtherNd(n)->pl().getGeom()) >
          6 * (g->getWidth(nf.edge) + g->getSpacing(nf.edge))) &&
        (nf.edge->getTo()->pl().getStops().size() == 0) &&
        (nf.edge->getFrom()->pl().getStops().size() == 0)) {
      res.push_back(nf);
    }
  }

  return res;
}

// _____________________________________________________________________________
std::set<NodeFront*> GraphBuilder::nodeGetOverlappingFronts(
    const RenderGraph* g, const LineNode* n) const {
  std::set<NodeFront*> ret;
  double minLength = 6;

  // TODO: why are nodefronts accessed via index?
  for (size_t i = 0; i < n->pl().getMainDirs().size(); ++i) {
    const NodeFront& fa = n->pl().getMainDirs()[i];

    for (size_t j = 0; j < n->pl().getMainDirs().size(); ++j) {
      const NodeFront& fb = n->pl().getMainDirs()[j];

      if (fa.geom.equals(fb.geom, 5) || j == i) continue;

      if (nodeFrontsOverlap(g, fa, fb)) {
        if (util::geo::len(*fa.edge->pl().getGeom()) > minLength &&
            fa.geom.distTo(*n->pl().getGeom()) <
                2.5 * g->getMaxNdFrontWidth(n)) {
          ret.insert(const_cast<NodeFront*>(&fa));
        }
        if (util::geo::len(*fb.edge->pl().getGeom()) > minLength &&
            fb.geom.distTo(*n->pl().getGeom()) <
                2.5 * g->getMaxNdFrontWidth(n)) {
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
  size_t numShr = g->getSharedRoutes(a.edge, b.edge).size();

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
      Ordering order(e->pl().getRoutes().size());
      for (size_t i = 0; i < e->pl().getRoutes().size(); i++) order[i] = i;
      c[e] = order;
    }
  }

  g->setConfig(c);
}
