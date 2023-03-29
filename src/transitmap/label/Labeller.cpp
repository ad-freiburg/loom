// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/label/Labeller.h"
#include "util/geo/Geo.h"

using shared::rendergraph::RenderGraph;
using transitmapper::label::Labeller;
using transitmapper::label::LineLabel;
using transitmapper::label::Overlaps;
using transitmapper::label::StationLabel;

using util::geo::MultiLine;
using util::geo::PolyLine;

// _____________________________________________________________________________
Labeller::Labeller(const config::Config* cfg) : _cfg(cfg) {}

// _____________________________________________________________________________
void Labeller::label(const RenderGraph& g, bool notDeg2) {
  labelStations(g, notDeg2);
  labelLines(g);
}

// _____________________________________________________________________________
util::geo::MultiLine<double> Labeller::getStationLblBand(
    const shared::linegraph::LineNode* n, double fontSize, uint8_t offset,
    const RenderGraph& g) {
  // TODO: the hull padding should be the same as in the renderer
  auto statHull = g.getStopGeoms(n, (_cfg->lineSpacing + _cfg->lineWidth) * 0.8,
                                 _cfg->tightStations, 4);

  double rad = util::geo::getEnclosingRadius(*n->pl().getGeom(), statHull);

  // TODO: determine the label width based on the real font width. This is
  // nontrivial, as it requires the fonts to be rendered for non-monospaced
  // fonts
  double labelW = (n->pl().stops().front().name.size() + 1) * fontSize / 2.1;

  util::geo::MultiLine<double> band;

  // TODO: should also be determined based on the font
  double h = fontSize * 0.75;

  util::geo::Line<double> geomBaseLine, geomMiddle, geomTop, capLeft, capRight;
  geomBaseLine.push_back(
      {n->pl().getGeom()->getX() + rad + (_cfg->lineSpacing + _cfg->lineWidth),
       n->pl().getGeom()->getY() - (offset * h / 2)});
  geomBaseLine.push_back({n->pl().getGeom()->getX() + rad + labelW,
                          n->pl().getGeom()->getY() - (offset * h / 2)});

  geomMiddle.push_back(
      {n->pl().getGeom()->getX() + rad + (_cfg->lineSpacing + _cfg->lineWidth),
       n->pl().getGeom()->getY() + h / 2 - (offset * h / 2)});
  geomMiddle.push_back({n->pl().getGeom()->getX() + rad + labelW,
                        n->pl().getGeom()->getY() + h / 2 - (offset * h / 2)});

  geomTop.push_back(
      {n->pl().getGeom()->getX() + rad + (_cfg->lineSpacing + _cfg->lineWidth),
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
void Labeller::labelStations(const RenderGraph& g, bool notdeg2) {
  std::vector<const shared::linegraph::LineNode*> orderedNds;
  for (auto n : g.getNds()) {
    if (n->pl().stops().size() == 0 || (notdeg2 && n->getDeg() == 2)) continue;
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

        auto overlaps = getOverlaps(band, n, g);

        if (overlaps.lineOverlaps + overlaps.statLabelOverlaps +
                overlaps.statOverlaps >
            0)
          continue;
        cands.push_back({band[0], band, fontSize, g.isTerminus(n), deg, offset,
                         overlaps, n->pl().stops().front()});
      }
    }

    std::sort(cands.begin(), cands.end());
    if (cands.size() == 0) continue;

    auto cand = cands.front();
    _stationLabels.push_back(cand);
    _statLblIdx.add(cand.band, _stationLabels.size() - 1);
  }
}

// _____________________________________________________________________________
Overlaps Labeller::getOverlaps(const util::geo::MultiLine<double>& band,
                               const shared::linegraph::LineNode* forNd,
                               const RenderGraph& g) const {
  std::set<const shared::linegraph::LineEdge*> proced;

  Overlaps ret{0, 0, 0, 0};

  std::set<const shared::linegraph::LineNode*> procedNds{forNd};

  for (auto line : band) {
    auto neighs = g.getNeighborEdges(
        line, g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing));
    for (auto neigh : neighs) {
      if (proced.count(neigh)) continue;

      if (util::geo::dist(*neigh->pl().getGeom(), band) <
          g.getTotalWidth(neigh) / 2) {
        ret.lineOverlaps++;
      }
      proced.insert(neigh);

      std::vector<const shared::linegraph::LineNode*> nds = {neigh->getTo(),
                                                             neigh->getFrom()};

      for (auto nd : nds) {
        if (nd->pl().stops().size() && !procedNds.count(nd)) {
          procedNds.insert(nd);

          auto statHull =
              g.getStopGeoms(nd, (_cfg->lineSpacing + _cfg->lineWidth) * 0.8,
                             _cfg->tightStations, 4);

          double rad =
              util::geo::getEnclosingRadius(*nd->pl().getGeom(), statHull);

          if (util::geo::dist(*nd->pl().getGeom(), band) <
              rad + (_cfg->lineWidth + _cfg->lineSpacing) / 2) {
            ret.statOverlaps++;
          }
        }
      }
    }
  }

  std::set<size_t> labelNeighs;
  _statLblIdx.get(band,
                   g.getMaxLineNum() * (_cfg->lineWidth + _cfg->lineSpacing),
                   &labelNeighs);

  for (auto id : labelNeighs) {
    auto labelNeigh = _stationLabels[id];
    if (util::geo::dist(labelNeigh.band, band) < 1) ret.statLabelOverlaps++;
  }

  return ret;
}

// _____________________________________________________________________________
void Labeller::labelLines(const RenderGraph& g) {
  LineLblIdx labelIdx = LineLblIdx();
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
          _statLblIdx.get(
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

          std::set<size_t> lineLabelNeighs;
          labelIdx.get(cand.getLine(),
                        20 * (_cfg->lineWidth + _cfg->lineSpacing),
                        &lineLabelNeighs);

          std::vector<const shared::linegraph::Line*> lines;
          for (auto lo : e->pl().getLines()) {
            lines.push_back(lo.line);
          }

          for (auto neighLabelId : lineLabelNeighs) {
            auto neighLabel = _lineLabels[neighLabelId];
            if (neighLabel.lines == lines &&
                util::geo::dist(cand.getLine(), neighLabel.geom.getLine()) <
                    20 * (_cfg->lineWidth + _cfg->lineSpacing)) {
              block = true;
              break;
            }
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
      labelIdx.add(cands.front().geom.getLine(), _lineLabels.size() - 1);
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
