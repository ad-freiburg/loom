// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_LABEL_LABELLER_H_
#define TRANSITMAP_LABEL_LABELLER_H_

#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "util/geo/Grid.h"

namespace transitmapper {
namespace label {

//const static std::vector<double> DEG_PENS = {0, 2.5, 10, 3, 0, 3, 10, 2};
// starting 90 deg
const static std::vector<double> DEG_PENS = {0, 2, 6, 3, 1, 3, 6, 2};

struct LineLabel {
  util::geo::PolyLine<double> geom;
  double centerDist;
  double fontSize;

  std::vector<const shared::linegraph::Line*> lines;
};

inline bool operator<(const LineLabel& a, const LineLabel& b) {
  return a.centerDist < b.centerDist;
}

struct Overlaps {
  size_t lineOverlaps;
  size_t lineLabelOverlaps;
  size_t statLabelOverlaps;
};

inline bool statNdCmp(const shared::linegraph::LineNode* a,
                      const shared::linegraph::LineNode* b) {
  // first degree 1 nodes
  size_t ad = a->getDeg();
  size_t bd = b->getDeg();
  if (ad == 1) ad = std::numeric_limits<size_t>::max();
  if (bd == 1) bd = std::numeric_limits<size_t>::max();
  return (ad > bd ||
          (ad == bd && shared::linegraph::LineGraph::getLDeg(a) >
                           shared::linegraph::LineGraph::getLDeg(b)));
}

struct StationLabel {
  util::geo::PolyLine<double> geom;
  util::geo::MultiLine<double> band;
  double fontSize;
  bool bold;

  size_t deg;
  size_t pos;
  Overlaps overlaps;

  shared::linegraph::Station s;

  double getPen() const {
    double score = overlaps.lineOverlaps * 15 +
                   overlaps.statLabelOverlaps * 20 +
                   overlaps.lineLabelOverlaps * 15;
    score += DEG_PENS[deg];

    if (pos == 1) score += 0.5;
    if (pos == 2) score += 0.1;
    return score;
  }
};

inline bool operator<(const StationLabel& a, const StationLabel& b) {
  return a.getPen() < b.getPen();
}

typedef util::geo::Grid<size_t, util::geo::MultiLine, double> StatLblGrid;

class Labeller {
 public:
  Labeller(const config::Config* cfg);

  void label(const shared::rendergraph::RenderGraph& g);

  const std::vector<LineLabel>& getLineLabels() const;
  const std::vector<StationLabel>& getStationLabels() const;

  util::geo::Box<double> getBBox() const;

 private:
  std::vector<LineLabel> _lineLabels;
  std::vector<StationLabel> _stationLabels;

  StatLblGrid _statLblGrid;

  const config::Config* _cfg;

  void labelStations(const shared::rendergraph::RenderGraph& g);
  void labelLines(const shared::rendergraph::RenderGraph& g);

  Overlaps getOverlaps(const util::geo::MultiLine<double>& band,
                       const shared::rendergraph::RenderGraph& g) const;

  util::geo::MultiLine<double> getStationLblBand(
      const shared::linegraph::LineNode* n, double fontSize, uint8_t offset,
      const shared::rendergraph::RenderGraph& g);
};
}  // namespace label
}  // namespace transitmapper

#endif  // TRANSITMAP_OUTPUT_SVGRENDERER_H_
