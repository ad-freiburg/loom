// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <istream>
#include <set>
#include <stack>
#include <vector>
#include "GraphBuilder.h"
#include "json/json.hpp"
#include "loom/config/TransitMapConfig.h"
#include "shared/rendergraph/RenderGraph.h"
#include "shared/linegraph/Line.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using loom::graph::GraphBuilder;
using shared::rendergraph::RenderGraph;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::Ordering;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;
using shared::linegraph::NodeFront;
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
