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
#include "topo/restr/RestrGraph.h"
#include "util/geo/Geo.h"
#include "util/geo/Grid.h"
#include "util/geo/PolyLine.h"
#include "util/graph/Graph.h"

using util::geo::Box;
using util::geo::DBox;
using util::geo::DLine;
using util::geo::DPoint;
using util::geo::Grid;
using util::geo::Line;
using util::geo::Point;
using util::geo::PolyLine;
using util::geo::SharedSegment;

using topo::config::TopoConfig;

using shared::linegraph::LineEdge;
using shared::linegraph::LineEdgePair;
using shared::linegraph::LineEdgePL;
using shared::linegraph::LineGraph;
using shared::linegraph::LineNode;
using shared::linegraph::LineNodePL;
using shared::linegraph::Station;

typedef Grid<LineNode*, Point, double> NodeGrid;

typedef std::map<const LineEdge*, std::set<const LineEdge*>> OrigEdgs;

namespace topo {

struct AggrDistFunc {
  AggrDistFunc(double maxDist) : _maxDist(maxDist){};

  double operator()(const LineEdge* a, const LineEdge* b) const {
    return fmax(_maxDist, a->pl().getLines().size() * (_maxDist / 2) +
                              b->pl().getLines().size() * (_maxDist / 2));
  };

  double operator()(const LineEdge* a) const { return operator()(a, a); };

  double _maxDist;
};

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

  int collapseShrdSegs();
  int collapseShrdSegs(double dCut);
  int collapseShrdSegs(double dCut, size_t MAX_ITERS);

  void averageNodePositions();
  void removeEdgeArtifacts();
  void removeNodeArtifacts(bool keepStations);

  void reconstructIntersections();

  size_t freeze();

  bool cleanUpGeoms();

  const OrigEdgs& freezeTrack(size_t i) const { return _origEdgs[i]; }
  void removeOrphanLines();

 private:
  const config::TopoConfig* _cfg;
  LineGraph* _g;

  LineNode* ndCollapseCand(const std::set<LineNode*>& notFrom,
                           const size_t numLines, const double maxD,
                           const util::geo::Point<double>& point,
                           const LineNode* spanA, const LineNode* spanB,
                           NodeGrid& grid, LineGraph* g) const;

  double maxD(size_t lines, const LineNode* nd, double d) const;
  double maxD(size_t lines, double d) const;
  double maxD(const LineNode* ndA, const LineNode* ndB, double d) const;

  bool combineNodes(LineNode* a, LineNode* b, LineGraph* g);
  bool combineEdges(LineEdge* a, LineEdge* b, LineNode* n, LineGraph* g);

  bool combineNodes(LineNode* a, LineNode* b);
  bool combineEdges(LineEdge* a, LineEdge* b, LineNode* n);

  void densifyEdg(LineEdge* e, LineGraph* g, double SEGL);

  bool contractNodes();

  void combContEdgs(const LineEdge* a, const LineEdge* b);
  void delOrigEdgsFor(const LineEdge* a);
  void delOrigEdgsFor(const LineNode* a);

  bool lineEq(const LineEdge* a, const LineEdge* b);

  bool contractEdges(bool keepStations);

  bool foldEdges(LineEdge* a, LineEdge* b);

  void mergeLines(LineEdge* newE, LineEdge* oldE, LineNode* newFrom,
                  LineNode* newTo);
  void supportEdge(LineEdge* ex, LineGraph* g);

  PolyLine<double> geomAvg(const LineEdgePL& geomA, double startA, double endA,
                           const LineEdgePL& geomB, double startB, double endB);

  DBox bbox() const;

  LineEdgePair split(LineEdgePL& a, LineNode* fr, LineNode* to, double p);

  std::set<const LineEdge*> _indEdges;
  std::set<LineEdgePair> _indEdgesPairs;
  std::map<LineEdgePair, size_t> _pEdges;

  std::vector<OrigEdgs> _origEdgs;
};

}  // namespace topo

#endif  // TOPO_MAPCONSTRUCTOR_MAPCONSTRUCTOR_H_
