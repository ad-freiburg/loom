// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_MAPCONSTRUCTOR_MAPCONSTRUCTOR_H_
#define TOPO_MAPCONSTRUCTOR_MAPCONSTRUCTOR_H_

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

using topo::config::TopoConfig;

using shared::transitgraph::TransitGraph;
using shared::transitgraph::TransitNode;
using shared::transitgraph::Station;
using shared::transitgraph::TransitEdge;
using shared::transitgraph::TransitEdgePair;
using shared::transitgraph::TransitNodePL;
using shared::transitgraph::TransitEdgePL;

typedef Grid<TransitNode*, Point, double> NodeGrid;
typedef Grid<TransitEdge*, Line, double> EdgeGrid;

typedef std::map<const TransitEdge*, std::set<const TransitEdge*>> OrigEdgs;

namespace topo {

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0){};
  ShrdSegWrap(TransitEdge* e, TransitEdge* f, SharedSegment<double> s)
      : e(e), f(f), s(s){};
  TransitEdge* e;
  TransitEdge* f;
  SharedSegment<double> s;
};

class MapConstructor {
 public:
  MapConstructor(const TopoConfig* cfg, TransitGraph* g);

  bool collapseShrdSegs();
  bool collapseShrdSegs(double dCut);
  bool collapseShrdSegs(double dCut, size_t steps);

  void averageNodePositions();
  void removeEdgeArtifacts();
  void removeNodeArtifacts();

  // TODO: make private
  bool contractNodes();

  size_t freeze();

  bool cleanUpGeoms();

  const OrigEdgs& freezeTrack(size_t i) const { return _origEdgs[i]; }

 private:
  const config::TopoConfig* _cfg;
  TransitGraph* _g;

  void routeDirRepl(TransitNode* oldN, TransitNode* newN, TransitEdge* e);

  ShrdSegWrap nextShrdSeg(double dCut, EdgeGrid* grid);
  bool combineNodes(TransitNode* a, TransitNode* b);
  bool combineEdges(TransitEdge* a, TransitEdge* b, TransitNode* n);

  void combContEdgs(const TransitEdge* a, const TransitEdge* b);

  bool routeEq(const TransitEdge* a, const TransitEdge* b);

  bool contractEdges();

  bool foldEdges(TransitEdge* a, TransitEdge* b);

  PolyLine<double> geomAvg(const TransitEdgePL& geomA, double startA,
                           double endA, const TransitEdgePL& geomB,
                           double startB, double endB);

  DBox bbox() const;
  EdgeGrid geoIndex();

  void edgeRpl(TransitNode* n, const TransitEdge* oldE,
               const TransitEdge* newE);

  TransitEdgePair split(TransitEdgePL& a, TransitNode* fr, TransitNode* to,
                        double p);

  std::set<const TransitEdge*> _indEdges;
  std::set<TransitEdgePair> _indEdgesPairs;
  std::map<TransitEdgePair, size_t> _pEdges;

  std::vector<OrigEdgs> _origEdgs;
};

}  // namespace topo

#endif  // TOPO_MAPCONSTRUCTOR_MAPCONSTRUCTOR_H_
