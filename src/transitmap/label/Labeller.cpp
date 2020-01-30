// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/label/Labeller.h"
#include "util/geo/Geo.h"

using transitmapper::label::Labeller;
using transitmapper::label::LineLabel;
using transitmapper::label::StationLabel;

using util::geo::PolyLine;

// _____________________________________________________________________________
void Labeller::label(const graph::RenderGraph& g) {
  labelStations(g);
  labelLines(g);
}

// _____________________________________________________________________________
void Labeller::labelStations(const graph::RenderGraph& g) {
  for (auto n : g.getNds()) {
    if (n->pl().stops().size() == 0) continue;

    util::geo::Line<double> geom;

    double fontSize =
        3 * (_cfg->lineWidth + _cfg->lineSpacing);

    double labelW = n->pl().stops().front().name.size() * fontSize / 1.5;

    geom.push_back(*n->pl().getGeom());
    geom.push_back({n->pl().getGeom()->getX() + labelW, n->pl().getGeom()->getY()});

    _stationLabels.push_back({geom, fontSize, n->pl().stops().front()});
  }
}

// _____________________________________________________________________________
void Labeller::labelLines(const graph::RenderGraph& g) {
  for (auto n : g.getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      double geomLen = util::geo::len(*e->pl().getGeom());

      // estimate label width
      double fontSize =
          3 * (_cfg->lineWidth + _cfg->lineSpacing);
      double labelW = ((fontSize / 3) * (e->pl().getLines().size() - 1));

      for (auto lo : e->pl().getLines()) {
        labelW += lo.line->label().size() * (fontSize);
      }

      // try out positions
      double step = geomLen / 4;

      std::vector<LineLabel> cands;

      for (int dir = -1; dir < 2; dir += 2) {
        double start = 0;
        while (start + labelW <= geomLen) {
          PolyLine<double> cand = util::geo::segment(
              *e->pl().getGeom(), start / geomLen, (start + labelW) / geomLen);
          cand.offsetPerp(dir * (g.getTotalWidth(e) / 2 +
                                 (_cfg->lineSpacing + _cfg->lineWidth)));

          bool block = false;

          for (auto neigh : g.getNeighborEdges(
                   cand.getLine(),
                   g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing) +
                       fontSize * 4)) {
            if (neigh == e) continue;
            if (util::geo::dist(cand.getLine(), *neigh->pl().getGeom()) <
                (g.getTotalWidth(neigh) / 2) + (fontSize)) {
              block = true;
              break;
            }
          }

          if (dir < 0) cand.reverse();

          std::vector<const shared::linegraph::Line*> lines;
          for (auto lo : e->pl().getLines()) {
            lines.push_back(lo.line);
          }

          if (!block)
            cands.push_back(
                {cand, fabs((geomLen / 2) - (start + (labelW / 2))), fontSize, lines});
          start += step;
        }
      }

      std::sort(cands.begin(), cands.end());
      if (cands.size() == 0) continue;
      _lineLabels.push_back(cands.front());
    }
  }
}

// _____________________________________________________________________________
const std::vector<LineLabel>& Labeller::getLineLabels() const {
  return _lineLabels;
}

// _____________________________________________________________________________
const std::vector<StationLabel>& Labeller::getStationLabels() const {
  return _stationLabels;
}
