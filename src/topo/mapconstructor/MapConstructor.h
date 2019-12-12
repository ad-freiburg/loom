// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_MAPCONSTRUCTOR_MAPCONSTRUCTOR_H_
#define TOPO_MAPCONSTRUCTOR_MAPCONSTRUCTOR_H_

#include <proj_api.h>
#include <algorithm>
#include <unordered_map>
#include "shared/linegraph/LineGraph.h"
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

using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::Station;
using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePair;
using shared::linegraph::LineNodePL;
using shared::linegraph::LineEdgePL;

typedef Grid<LineNode*, Point, double> NodeGrid;
typedef Grid<LineEdge*, Line, double> EdgeGrid;

typedef std::map<const LineEdge*, std::set<const LineEdge*>> OrigEdgs;

namespace topo {

struct ShrdSegWrap {
  ShrdSegWrap() : e(0), f(0){};
  ShrdSegWrap(LineEdge* e, LineEdge* f, SharedSegment<double> s)
      : e(e), f(f), s(s){};
  LineEdge* e;
  LineEdge* f;
  SharedSegment<double> s;
};

class MapConstructor {
 public:
  MapConstructor(const TopoConfig* cfg, LineGraph* g);

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
  LineGraph* _g;

  void routeDirRepl(LineNode* oldN, LineNode* newN, LineEdge* e);

  ShrdSegWrap nextShrdSeg(double dCut, EdgeGrid* grid);
  bool combineNodes(LineNode* a, LineNode* b);
  bool combineEdges(LineEdge* a, LineEdge* b, LineNode* n);

  void combContEdgs(const LineEdge* a, const LineEdge* b);

  bool routeEq(const LineEdge* a, const LineEdge* b);

  bool contractEdges();

  bool foldEdges(LineEdge* a, LineEdge* b);

  PolyLine<double> geomAvg(const LineEdgePL& geomA, double startA, double endA,
                           const LineEdgePL& geomB, double startB, double endB);

  DBox bbox() const;
  EdgeGrid geoIndex();

  void edgeRpl(LineNode* n, const LineEdge* oldE, const LineEdge* newE);

  LineEdgePair split(LineEdgePL& a, LineNode* fr, LineNode* to, double p);

  std::set<const LineEdge*> _indEdges;
  std::set<LineEdgePair> _indEdgesPairs;
  std::map<LineEdgePair, size_t> _pEdges;

  std::vector<OrigEdgs> _origEdgs;
};

}  // namespace topo

#endif  // TOPO_MAPCONSTRUCTOR_MAPCONSTRUCTOR_H_
