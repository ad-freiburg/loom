// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_BUILDER_BUILDER_H_
#define TOPO_BUILDER_BUILDER_H_

#include <proj_api.h>
#include <algorithm>
#include <unordered_map>
#include "shared/transitgraph/TransitGraph.h"
#include "topo/config/TopoConfig.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Graph.h"

using util::geo::Grid;
using util::geo::Box;
using util::geo::Line;
using util::geo::PolyLine;
using util::geo::Point;
using util::geo::Line;
using util::geo::DPoint;
using util::geo::DLine;
using util::geo::DBox;
using util::geo::SharedSegment;

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::TransitEdge;
using shared::transitgraph::TransitNodePL;
using shared::transitgraph::TransitEdgePL;

typedef Grid<TransitNode*, Point, double> NodeGrid;
typedef Grid<TransitEdge*, Line, double> EdgeGrid;

namespace topo {

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0){};
  ShrdSegWrap(TransitEdge* e, TransitEdge* f, SharedSegment<double> s)
      : e(e), f(f), s(s){};
  TransitEdge* e;
  TransitEdge* f;
  SharedSegment<double> s;
};

class Builder {
 public:
  Builder(const config::TopoConfig* cfg);

  bool createTopologicalNodes(TransitGraph* g, bool final);
  bool createTopologicalNodes(TransitGraph* g, bool final, size_t steps);
  void averageNodePositions(TransitGraph* g);
  void removeEdgeArtifacts(TransitGraph* g);
  void removeNodeArtifacts(TransitGraph* g);
  void cleanEx(TransitGraph* tg) const;

  // TODO: make private
  bool contractNodes(TransitGraph* g);

  void inferRestrictions(TransitGraph* g) const;

 private:
  const config::TopoConfig* _cfg;
  projPJ _mercProj;
  projPJ _graphProj;

  void routeDirRepl(TransitNode* oldN, TransitNode* newN, TransitEdge* e) const;

  ShrdSegWrap getNextSharedSegment(TransitGraph* g, bool final,
                                   EdgeGrid* grid) const;
  PolyLine<double> getAveragedFromSharedSeg(const ShrdSegWrap& w) const;

  bool combineNodes(TransitNode* a, TransitNode* b, TransitGraph* g);
  bool combineEdges(TransitEdge* a, TransitEdge* b, TransitNode* n,
                    TransitGraph* g);

  void terminusPass(TransitNode* b, const TransitEdge* connecting);

  bool crossesAt(const TransitNode* a, const TransitEdge* e,
                 const TransitEdge* f) const;

  bool routeEq(const TransitEdge* a, const TransitEdge* b) const;

  bool contractEdges(TransitGraph* g);

  bool foldEdges(TransitEdge* a, TransitEdge* b);

  bool isTriFace(const TransitEdge* a, const TransitGraph* g) const;

  PolyLine<double> geomAvg(const TransitEdgePL& geomA, double startA,
                           double endA, const TransitEdgePL& geomB,
                           double startB, double endB) const;

  DBox getGraphBoundingBox(const TransitGraph* g) const;
  EdgeGrid getGeoIndex(const TransitGraph* g) const;

  std::pair<TransitEdge*, TransitEdge*> split(TransitEdgePL& a, TransitNode* fr,
                                              TransitNode* to,
                                              TransitGraph* g) const;

  mutable std::set<const TransitEdge*> _indEdges;
  mutable std::set<std::pair<const TransitEdge*, const TransitEdge*> >
      _indEdgesPairs;
  mutable std::map<std::pair<const TransitEdge*, const TransitEdge*>, size_t>
      _pEdges;
};

}  // namespace topo

#endif  // TOPO_BUILDER_BUILDER_H_
