// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "transitmap/label/Labeller.h"
#include "util/geo/Geo.h"

using transitmapper::label::Labeller;
using transitmapper::label::LineLabel;
using transitmapper::label::StationLabel;
using transitmapper::label::Overlaps;

using util::geo::PolyLine;
using util::geo::MultiLine;

// _____________________________________________________________________________
Labeller::Labeller(const config::Config* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void Labeller::label(const graph::RenderGraph& g) {
  // leave enough room for labels
  auto bbox = util::geo::pad(g.getBBox(), 500);
  _statLblGrid = StatLblGrid(200, 200, bbox);

  labelStations(g);
  labelLines(g);
}

// _____________________________________________________________________________
util::geo::MultiLine<double> Labeller::getStationLblBand(
    const shared::linegraph::LineNode* n, double fontSize, uint8_t offset,
    const graph::RenderGraph& g) {
  // TODO: the hull padding should be the same as in the renderer
  auto statHull =
      g.getStopGeoms(n, (_cfg->lineSpacing + _cfg->lineWidth) * 0.8,
                       _cfg->simpleRenderForTwoEdgeNodes);

  double rad = util::geo::getEnclosingRadius(*n->pl().getGeom(), statHull);

  double labelW = n->pl().stops().front().name.size() * fontSize / 1.6;

  util::geo::MultiLine<double> band;

  double h = fontSize * 0.75;

  util::geo::Line<double> geomBaseLine, geomMiddle, geomTop, capLeft, capRight;
  geomBaseLine.push_back({n->pl().getGeom()->getX() + rad +
                              (_cfg->lineSpacing + _cfg->lineWidth),
                          n->pl().getGeom()->getY() - (offset * h / 2)});
  geomBaseLine.push_back({n->pl().getGeom()->getX() + rad + labelW,
                          n->pl().getGeom()->getY() - (offset * h / 2)});

  geomMiddle.push_back({n->pl().getGeom()->getX() + rad +
                            (_cfg->lineSpacing + _cfg->lineWidth),
                        n->pl().getGeom()->getY() + h / 2 - (offset * h / 2)});
  geomMiddle.push_back({n->pl().getGeom()->getX() + rad + labelW,
                        n->pl().getGeom()->getY() + h / 2 - (offset * h / 2)});

  geomTop.push_back({n->pl().getGeom()->getX() + rad +
                         (_cfg->lineSpacing + _cfg->lineWidth),
                     n->pl().getGeom()->getY() + h - (offset * h / 2)});
  geomTop.push_back({n->pl().getGeom()->getX() + rad + labelW,
                     n->pl().getGeom()->getY() + h - (offset * h / 2)});

  capLeft = util::geo::PolyLine<double>(geomMiddle)
                .getOrthoLineAtDist(util::geo::len(geomMiddle), h)
                .getLine();
  capRight = util::geo::PolyLine<double>(geomMiddle)
                 .getOrthoLineAtDist(0, h)
                 .getLine();

  band.push_back(geomBaseLine);
  band.push_back(geomMiddle);
  band.push_back(geomTop);
  band.push_back(capLeft);
  band.push_back(capRight);

  return band;
}

// _____________________________________________________________________________
void Labeller::labelStations(const graph::RenderGraph& g) {
  std::vector<const shared::linegraph::LineNode*> orderedNds;
  for (auto n : g.getNds()) {
    if (n->pl().stops().size() == 0) continue;
    orderedNds.push_back(n);
  }

  std::sort(orderedNds.begin(), orderedNds.end(), statNdCmp);

  for (auto n : orderedNds) {

    double fontSize = _cfg->stationLabelSize;

    std::vector<StationLabel> cands;

    for (uint8_t offset = 0; offset < 3; offset++) {
      for (size_t deg = 0; deg < 8; deg++) {
        auto band = getStationLblBand(n, fontSize, offset, g);
        band = util::geo::rotate(band, 45 * deg, *n->pl().getGeom());

        auto overlaps = getOverlaps(band, g);

        if (overlaps.lineOverlaps + overlaps.statLabelOverlaps > 0) continue;
        cands.push_back({band[0], band, fontSize, g.isTerminus(n), deg, offset,
                         overlaps, n->pl().stops().front()});
      }
    }

    std::sort(cands.begin(), cands.end());
    if (cands.size() == 0) continue;

    auto cand = cands.front();
    _stationLabels.push_back(cand);
    _statLblGrid.add(cand.band, _stationLabels.size() - 1);
  }
}

// _____________________________________________________________________________
Overlaps Labeller::getOverlaps(const util::geo::MultiLine<double>& band,
                               const graph::RenderGraph& g) const {
  std::set<const shared::linegraph::LineEdge*> proced;

  Overlaps ret{0, 0, 0};

  for (auto line : band) {
    for (auto neigh : g.getNeighborEdges(
             line, g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing))) {
      if (proced.count(neigh)) continue;

      if (util::geo::dist(*neigh->pl().getGeom(), band) <
          g.getTotalWidth(neigh) / 2) {
        ret.lineOverlaps++;
      }
      proced.insert(neigh);
    }
  }

  std::set<size_t> labelNeighs;
  _statLblGrid.get(band,
                   g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing),
                   &labelNeighs);

  for (auto id : labelNeighs) {
    auto labelNeigh = _stationLabels[id];
    if (util::geo::dist(labelNeigh.band, band) < 1) {
      ret.statLabelOverlaps++;
    }
  }

  return ret;
}

// _____________________________________________________________________________
void Labeller::labelLines(const graph::RenderGraph& g) {
  for (auto n : g.getNds()) {
    for (auto e : n->getAdjList()) {
      if (e->getFrom() != n) continue;
      double geomLen = util::geo::len(*e->pl().getGeom());

      // estimate label width
      double fontSize = _cfg->lineLabelSize;

      double labelW = ((fontSize / 3) * (e->pl().getLines().size() - 1));

      for (auto lo : e->pl().getLines()) {
        labelW += lo.line->label().size() * (fontSize);
      }

      // try out positions
      double step = fontSize;

      std::vector<LineLabel> cands;

      for (int dir = -1; dir < 2; dir += 2) {
        double start = 0;
        while (start + labelW <= geomLen) {
          PolyLine<double> cand = util::geo::segment(
              *e->pl().getGeom(), start / geomLen, (start + labelW) / geomLen);
          if (cand.getLength() < 5) break;
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

          std::set<size_t> labelNeighs;
          _statLblGrid.get(
              MultiLine<double>{cand.getLine()},
              g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing),
              &labelNeighs);

          for (auto neighId : labelNeighs) {
            auto neigh = _stationLabels[neighId];
            if (util::geo::dist(cand.getLine(), neigh.band) < (fontSize)) {
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
            cands.push_back({cand, fabs((geomLen / 2) - (start + (labelW / 2))),
                             fontSize, lines});
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

// _____________________________________________________________________________
util::geo::Box<double> Labeller::getBBox() const {
  util::geo::Box<double> ret;

  for (auto lbl : _lineLabels)
    ret = util::geo::extendBox(lbl.geom.getLine(), ret);
  for (auto lbl : _stationLabels) ret = util::geo::extendBox(lbl.band, ret);

  return ret;
}
