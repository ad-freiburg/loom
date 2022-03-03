// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TOPO_STATINSERTER_STATINSERTER_H_
#define TOPO_STATINSERTER_STATINSERTER_H_

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

struct StationOcc {
  Station station;
  std::set<const LineEdge*> edges;
};

struct StationCand {
  // either an edge
  LineEdge* edg;
  double pos;

  // or a node
  LineNode* nd;

  double dist;

  size_t shouldServ;

  size_t truelyServ;
  size_t falselyServ;
};

class StatInserter {
 public:
  StatInserter(const TopoConfig* cfg, LineGraph* g);

  void init();
  bool insertStations(const OrigEdgs& origEdgs);

 private:
  const config::TopoConfig* _cfg;
  LineGraph* _g;

  std::vector<StationCand> candidates(const StationOcc& occ,
                                      const EdgeGrid& idx,
                                      const OrigEdgs& origEdgs);

  DBox bbox() const;
  EdgeGrid geoIndex();

  static double candScore(const StationCand& c);

  std::pair<size_t, size_t> served(const std::vector<LineEdge*>& adj,
                                   const std::set<const LineEdge*>& toServe,
                                   const OrigEdgs& origEdgs);

  LineEdgePair split(LineEdgePL& a, LineNode* fr, LineNode* to, double p);

  void edgeRpl(LineNode* n, const LineEdge* oldE, const LineEdge* newE);

  std::vector<std::vector<StationOcc>> _statClusters;
};

}  // namespace topo

#endif  // TOPO_STATINSERTER_STATINSERTER_H_
