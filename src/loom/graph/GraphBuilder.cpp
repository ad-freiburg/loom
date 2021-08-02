// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <istream>
#include <set>
#include <stack>
#include <vector>
#include "GraphBuilder.h"
#include "json/json.hpp"
#include "loom/config/LoomConfig.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using loom::graph::GraphBuilder;
using shared::linegraph::Line;
using shared::linegraph::LineEdge;
using shared::linegraph::LineNode;
using shared::linegraph::NodeFront;
using shared::rendergraph::OrderCfg;
using shared::rendergraph::Ordering;
using shared::rendergraph::RenderGraph;
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

      n->pl().addFront(f);
    }
  }
}
