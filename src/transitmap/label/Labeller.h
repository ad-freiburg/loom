// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_LABEL_LABELLER_H_
#define TRANSITMAP_LABEL_LABELLER_H_

#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/RenderGraph.h"
#include "shared/linegraph/Line.h"

namespace transitmapper {
namespace label {

struct LineLabel {
  util::geo::PolyLine<double> geom;
  double centerDist;
  double fontSize;

  std::vector<const shared::linegraph::Line*> lines;
};

inline bool operator<(const LineLabel& a, const LineLabel& b) {
  return a.centerDist < b.centerDist;
}

struct StationLabel {
  util::geo::PolyLine<double> geom;
  double fontSize;

  shared::linegraph::Station s;
};

class Labeller {
 public:
  Labeller(const config::Config* cfg) : _cfg(cfg) {}

  void label(const graph::RenderGraph& g);

  const std::vector<LineLabel>& getLineLabels() const;
  const std::vector<StationLabel>& getStationLabels() const;

 private:
  std::vector<LineLabel> _lineLabels;
  std::vector<StationLabel> _stationLabels;

  const config::Config* _cfg;

  void labelStations(const graph::RenderGraph& g);
  void labelLines(const graph::RenderGraph& g);
};

}
}

#endif  // TRANSITMAP_OUTPUT_SVGRENDERER_H_
