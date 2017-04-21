// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>
#include <ostream>
#include "./../config/TransitMapConfig.h"
#include "./../graph/Route.h"
#include "./SvgOutput.h"
#include "pbutil/String.h"
#include "pbutil/geo/PolyLine.h"

using namespace transitmapper;
using namespace output;

using pbutil::toString;

// _____________________________________________________________________________
SvgOutput::SvgOutput(std::ostream* o, const config::Config* cfg)
    : _o(o), _w(o, true), _cfg(cfg) {}

// _____________________________________________________________________________
void SvgOutput::print(const graph::TransitGraph& outG) {
  std::map<std::string, std::string> params;
  RenderParams rparams;

  rparams.xOff = outG.getBoundingBox().min_corner().get<0>();
  rparams.yOff = outG.getBoundingBox().min_corner().get<1>();

  rparams.width = outG.getBoundingBox().max_corner().get<0>() - rparams.xOff;
  rparams.height = outG.getBoundingBox().max_corner().get<1>() - rparams.yOff;

  rparams.width *= _cfg->outputResolution;
  rparams.height *= _cfg->outputResolution;

  params["width"] = std::to_string(rparams.width) + "px";
  params["height"] = std::to_string(rparams.height) + "px";

  *_o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  *_o << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
         "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">";

  outputEdges(outG, rparams);

  _w.openTag("svg", params);

  _w.openTag("defs");

  for (auto const& m : _markers) {
    params.clear();
    params["id"] = m.name;
    params["orient"] = "auto";
    params["markerWidth"] = "2";
    params["markerHeight"] = "1";
    params["refY"] = "0.5";
    params["refX"] = "0";

    _w.openTag("marker", params);

    params.clear();
    params["d"] = m.path;
    params["fill"] = m.color;
    params["stroke"] = "";

    _w.openTag("path", params);

    _w.closeTag();
    _w.closeTag();

    params.clear();
    params["id"] = m.name + "_black";
    params["orient"] = "auto";
    params["markerWidth"] = "2";
    params["markerHeight"] = "1";
    params["refY"] = "0.5";
    params["refX"] = "0";

    _w.openTag("marker", params);

    params.clear();
    params["d"] = m.path;
    params["fill"] = "#000000";

    _w.openTag("path", params);

    _w.closeTag();
    _w.closeTag();
  }

  _w.closeTag();

  for (graph::Node* n : outG.getNodes()) {
    renderNodeConnections(outG, n, rparams);
  }

  renderDelegates(outG, rparams);

  for (graph::Node* n : outG.getNodes()) {
    if (_cfg->renderStationNames) renderNodeScore(outG, n, rparams);
  }

  outputNodes(outG, rparams);
  if (_cfg->renderNodeFronts) {
    renderNodeFronts(outG, rparams);
  }

  _w.closeTags();
}

// _____________________________________________________________________________
void SvgOutput::outputNodes(const graph::TransitGraph& outG,
                            const RenderParams& rparams) {
  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    std::map<std::string, std::string> params;

    if (_cfg->renderStations && n->getStops().size() > 0 &&
        n->getMainDirs().size() > 0) {
      params["stroke"] = "black";
      params["stroke-width"] = "1";
      params["fill"] = "white";
      printPolygon(n->getConvexFrontHull(20), params, rparams);
    } else if (false) {
      params["r"] = "5";
      params["fill"] = "#FF00FF";
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::renderNodeFronts(const graph::TransitGraph& outG,
                                 const RenderParams& rparams) {
  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    for (auto& f : n->getMainDirs()) {
      const PolyLine p = f.geom;
      std::stringstream style;
      style << "fill:none;stroke:red"
            << ";stroke-linejoin: "
               "miter;stroke-linecap:round;stroke-opacity:0.5;stroke-width:1";
      std::map<std::string, std::string> params;
      params["style"] = style.str();
      printLine(p, params, rparams);

      Point a = p.getPointAt(.5).p;

      std::stringstream styleA;
      styleA << "fill:none;stroke:red"
             << ";stroke-linejoin: "
                "miter;stroke-linecap:round;stroke-opacity:1;stroke-width:.5";
      params["style"] = styleA.str();

      printLine(PolyLine(n->getPos(), a), params, rparams);
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::outputEdges(const graph::TransitGraph& outG,
                            const RenderParams& rparams) {
  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      renderEdgeTripGeom(outG, e, rparams);
    }
  }
}

// _____________________________________________________________________________
void SvgOutput::renderNodeConnections(const graph::TransitGraph& outG,
                                      const graph::Node* n,
                                      const RenderParams& rparams) {
  if (n->getStops().size() != 0) return;

  for (auto& ie :
       n->getInnerGeometries(outG.getConfig(), _cfg->innerGeometryPrecision)) {
    renderLinePart(ie.geom, ie.from.edge->getWidth(), *ie.from.route, 0);
  }
}

// _____________________________________________________________________________
void SvgOutput::renderLinePart(const PolyLine p, double width,
                               const graph::Route& route,
                               const Nullable<style::LineStyle> style) {
  renderLinePart(p, width, route, "", style);
}

// _____________________________________________________________________________
void SvgOutput::renderLinePart(const PolyLine p, double width,
                               const graph::Route& route,
                               const std::string& endMarker,
                               const Nullable<style::LineStyle> style) {
  std::stringstream styleOutline;
  styleOutline << "fill:none;stroke:#000000";

  if (!endMarker.empty()) {
    styleOutline << ";marker-end:url(#" << endMarker << "_black)";
  }

  if (!style.isNull()) {
    if (style.get().getDashArray().size()) {
      styleOutline << ";stroke-dasharray:" << style.get().getDashArrayString();
    }
  }

  styleOutline << ";stroke-linecap:round;stroke-opacity:0.8;stroke-width:"
               << (width + _cfg->outlineWidth) * _cfg->outputResolution;
  Params paramsOutline;
  paramsOutline["style"] = styleOutline.str();

  std::stringstream styleStr;
  styleStr << "fill:none;stroke:#" << route.color;

  if (!endMarker.empty()) {
    styleStr << ";marker-end:url(#" << endMarker << ")";
  }

  if (!style.isNull()) {
    if (style.get().getDashArray().size()) {
      styleStr << ";stroke-dasharray:" << style.get().getDashArrayString();
    }

    if (!style.get().getCss().empty()) {
      std::string css = style.get().getCss();
      pbutil::replaceAll(css, "\"", "&quot;");
      styleStr << ";" << css << ";";
    }
  }

  styleStr << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
           << width * _cfg->outputResolution;
  Params params;
  params["style"] = styleStr.str();

  _delegates[(uintptr_t)&route].push_back(OutlinePrintPair(
      PrintDelegate(paramsOutline, p), PrintDelegate(params, p)));
}

// _____________________________________________________________________________
void SvgOutput::renderNodeScore(const graph::TransitGraph& outG,
                                const graph::Node* n,
                                const RenderParams& rparams) {
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  Params params;
  params["x"] = std::to_string(
      (n->getPos().get<0>() - xOffset) * _cfg->outputResolution + 0);
  params["y"] = std::to_string(
      rparams.height -
      (n->getPos().get<1>() - yOffset) * _cfg->outputResolution - 0);
  params["style"] =
      "font-family:Verdana;font-size:8px; font-style:normal; font-weight: "
      "normal; fill: white; stroke-width: 0.25px; stroke-linecap: butt; "
      "stroke-linejoin: miter; stroke: black";
  _w.openTag("text", params);
  if (n->getStops().size()) {
    _w.writeText((*n->getStops().begin()).id);
    _w.writeText("\n");
  }

  _w.writeText(pbutil::toString(n));
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::renderEdgeTripGeom(const graph::TransitGraph& outG,
                                   const graph::Edge* e,
                                   const RenderParams& rparams) {
  const graph::NodeFront* nfTo = e->getTo()->getNodeFrontFor(e);
  const graph::NodeFront* nfFrom = e->getFrom()->getNodeFrontFor(e);

  PolyLine center = e->getGeom();

  double lineW = e->getWidth();
  double lineSpc = e->getSpacing();
  double offsetStep = lineW + lineSpc;
  double oo = e->getTotalWidth();

  double o = oo;

  assert(outG.getConfig().find(e) != outG.getConfig().end());

  size_t a = 0;
  for (size_t i : outG.getConfig().find(e)->second) {
    const graph::RouteOccurance& r = e->getTripsUnordered()[i];
    PolyLine p = center;

    double offset = -(o - oo / 2.0 - e->getWidth() / 2.0);

    p.offsetPerp(offset);

    std::set<LinePoint, LinePointCmp> iSects = nfTo->geom.getIntersections(p);
    if (iSects.size() > 0) {
      p = p.getSegment(0, iSects.begin()->totalPos);
    } else {
      p << nfTo->geom.projectOn(p.getLine().back()).p;
    }

    std::set<LinePoint, LinePointCmp> iSects2 =
        nfFrom->geom.getIntersections(p);
    if (iSects2.size() > 0) {
      p = p.getSegment(iSects2.begin()->totalPos, 1);
    } else {
      p >> nfFrom->geom.projectOn(p.getLine().front()).p;
    }

    double arrowLength = (6 / _cfg->outputResolution);

    if (r.direction != 0 && center.getLength() > arrowLength * 4) {
      std::stringstream markerName;
      markerName << e << ":" << r.route << ":" << i;

      std::string markerPathMale = getMarkerPathMale(lineW);
      std::string markerPathFemale = getMarkerPathFemale(lineW);
      EndMarker emm(markerName.str() + "_m", "#" + r.route->color,
                    markerPathMale, lineW, lineW);
      EndMarker emf(markerName.str() + "_f", "#" + r.route->color,
                    markerPathFemale, lineW, lineW);
      _markers.push_back(emm);
      _markers.push_back(emf);
      PolyLine firstPart =
          p.getSegmentAtDist(0, p.getLength() / 2 - arrowLength / 2);
      PolyLine secondPart = p.getSegmentAtDist(
          p.getLength() / 2 + arrowLength / 2, p.getLength());

      if (r.direction == e->getTo()) {
        renderLinePart(firstPart, lineW, *r.route, markerName.str() + "_m",
                       r.style);
        renderLinePart(secondPart.getReversed(), lineW, *r.route,
                       markerName.str() + "_f", r.style);
      } else {
        renderLinePart(firstPart, lineW, *r.route, markerName.str() + "_f",
                       r.style);
        renderLinePart(secondPart.getReversed(), lineW, *r.route,
                       markerName.str() + "_m", r.style);
      }
    } else {
      renderLinePart(p, lineW, *r.route, r.style);
    }

    a++;

    // break;
    o -= offsetStep;
  }
}

// _____________________________________________________________________________
std::string SvgOutput::getMarkerPathMale(double w) const {
  std::stringstream path;
  path << "M0,0 V1 H.5 L1.3,.5 L.5,0 Z";

  return path.str();
}

// _____________________________________________________________________________
std::string SvgOutput::getMarkerPathFemale(double w) const {
  std::stringstream path;
  path << "M1.3,1 L.5,.5 L1.3,0 L0,0 L0,1";

  return path.str();
}

// _____________________________________________________________________________
void SvgOutput::renderDelegates(const graph::TransitGraph& outG,
                                const RenderParams& rparams) {
  for (auto& a : _delegates) {
    _w.openTag("g");
    for (auto& pd : a.second) {
      printLine(pd.first.second, pd.first.first, rparams);
    }
    _w.closeTag();

    _w.openTag("g");
    for (auto& pd : a.second) {
      printLine(pd.second.second, pd.second.first, rparams);
    }
    _w.closeTag();
  }
}

// _____________________________________________________________________________
void SvgOutput::printPoint(const Point& p, const std::string& style,
                           const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["cx"] =
      std::to_string((p.get<0>() - rparams.xOff) * _cfg->outputResolution);
  params["cy"] = std::to_string(
      rparams.height - (p.get<1>() - rparams.yOff) * _cfg->outputResolution);
  params["r"] = "2";
  params["style"] = style;
  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::printLine(const PolyLine& l, const std::string& style,
                          const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printLine(l, params, rparams);
}

// _____________________________________________________________________________
void SvgOutput::printLine(const PolyLine& l,
                          const std::map<std::string, std::string>& ps,
                          const RenderParams& rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  for (auto& p : l.getLine()) {
    points << " " << (p.get<0>() - rparams.xOff) * _cfg->outputResolution << ","
           << rparams.height -
                  (p.get<1>() - rparams.yOff) * _cfg->outputResolution;
  }

  params["points"] = points.str();

  _w.openTag("polyline", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::printPolygon(const Polygon& g,
                             const std::map<std::string, std::string>& ps,
                             const RenderParams& rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  for (auto& p : g.outer()) {
    points << " " << (p.get<0>() - rparams.xOff) * _cfg->outputResolution << ","
           << rparams.height -
                  (p.get<1>() - rparams.yOff) * _cfg->outputResolution;
  }

  params["points"] = points.str();

  _w.openTag("polygon", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::printPolygon(const Polygon& g, const std::string& style,
                             const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printPolygon(g, params, rparams);
}
