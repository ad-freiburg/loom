// Copy//right 2016, U  //niversity of Freiburg,
//{/ Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>
#include <ostream>
#include "./../config/TransitMapConfig.h"
#include "./SvgOutput.h"
#include "../geo/PolyLine.h"

using namespace transitmapper;
using namespace output;


// _____________________________________________________________________________
SvgOutput::SvgOutput(std::ostream* o, const config::Config* cfg)
: _o(o), _w(o, true), _cfg(cfg) {

}

// _____________________________________________________________________________
void SvgOutput::print(const graph::TransitGraph& outG) {
  std::map<std::string, std::string> params;

  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();


  int64_t width = outG.getBoundingBox().max_corner().get<0>() - xOffset;
  int64_t height = outG.getBoundingBox().max_corner().get<1>() - yOffset;

  width *= _cfg->outputResolution;
  height *= _cfg->outputResolution;

  params["width"] = std::to_string(width) + "px";
  params["height"] = std::to_string(height) + "px";

  *_o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  *_o << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">";

  _w.openTag("svg", params);

  // TODO: output edges

  outputEdges(outG, width, height);

  for (graph::Node* n : outG.getNodes()) {
    renderNodeConnections(outG, n, width, height);
  }

  renderDelegates(outG, width, height);

  for (graph::Node* n : outG.getNodes()) {
    if (_cfg->renderStationNames) renderNodeScore(outG, n, width, height);
  }

  //outputNodes(outG, width, height);
  if (_cfg->renderNodeFronts) {
    renderNodeFronts(outG, width, height);
  }

  _w.closeTags();
}

// _____________________________________________________________________________
void SvgOutput::outputNodes(const graph::TransitGraph& outG, double w, double h)
{
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();


  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    std::map<std::string, std::string> params;

    if (_cfg->renderStations && n->getStops().size() > 0 && n->getMainDirs().size() > 0) {
      params["stroke"] = "black";
      params["stroke-width"] = "1";
      params["fill"] = "white";
      printPolygon(n->getConvexFrontHull(20), params, w, h, xOffset, yOffset);
    } else if (false) {
      params["r"] = "5";
      params["fill"] = "#FF00FF";
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::renderNodeFronts(const graph::TransitGraph& outG, double w,
    double h) {
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    for (auto& f : n->getMainDirs()) {
      const geo::PolyLine p = f.geom;
      std::stringstream style;
      style << "fill:none;stroke:red"
        << ";stroke-linecap:round;stroke-opacity:0.5;stroke-width:1";
      std::map<std::string, std::string> params;
      params["style"] = style.str();
      printLine(p, params, w, h, xOffset, yOffset);

      util::geo::Point dir = p.getPointAt(1).p;
      style.clear();
      style << "fill:none;stroke:red"
        << ";stroke-linecap:round;stroke-opacity:1;stroke-width:.5";
      printPoint(dir, style.str(), w, h, xOffset, yOffset);

      util::geo::Point a = p.getPointAt(.5).p;

      std::stringstream styleA;
      styleA << "fill:none;stroke:red"
        << ";stroke-linecap:round;stroke-opacity:1;stroke-width:.5";
      params["style"] = styleA.str();

      printLine(geo::PolyLine(n->getPos(), a), params, w, h, xOffset, yOffset);
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::outputEdges(const graph::TransitGraph& outG, double w, double h)
{

  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (const graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        renderEdgeTripGeom(outG, g, e, w, h);
      }
    }
  }
}

// _____________________________________________________________________________
void SvgOutput::renderNodeConnections(const graph::TransitGraph& outG,
    const graph::Node* n, double w, double h) {
  if (n->getStops().size() != 0) return;

  for (auto& ie : n->getInnerGeometries(outG.getConfig(), true)) {
    renderLinePart(ie.geom, ie.etg->getWidth(), *ie.route);
  }
}

// _____________________________________________________________________________
void SvgOutput::renderLinePart(const geo::PolyLine p, double width,
    const gtfs::Route& route) {
  std::stringstream styleOutline;
  styleOutline << "fill:none;stroke:#000000"
    << ";stroke-linecap:round;stroke-opacity:0.8;stroke-width:"
    << (width + 4) * _cfg->outputResolution;
  Params paramsOutline;
  paramsOutline["style"] = styleOutline.str();

  std::stringstream style;
  style << "fill:none;stroke:#" << route.getColorString()
    << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
    << width * _cfg->outputResolution;
  Params params;
  params["style"] = style.str();

  _delegates[(uintptr_t)&route].push_back(
      OutlinePrintPair(
        PrintDelegate(paramsOutline, p),
        PrintDelegate(params, p)));
}

// _____________________________________________________________________________
void SvgOutput::renderNodeScore(const graph::TransitGraph& outG,
    const graph::Node* n, double w, double h) {
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  Params params;
  params["x"] = std::to_string((n->getPos().get<0>() - xOffset) * _cfg->outputResolution + 0);
  params["y"] = std::to_string(h-(n->getPos().get<1>() - yOffset) * _cfg->outputResolution - 0);
  params["style"] = "font-family:Verdana;font-size:8px; font-style:normal; font-weight: normal; fill: white; stroke-width: 0.25px; stroke-linecap: butt; stroke-linejoin: miter; stroke: black";
  _w.openTag("text", params);
  if (n->getStops().size()) {
    _w.writeText((*n->getStops().begin())->getId());
    _w.writeText("\n");
  }
  //_w.writeText(std::to_string(n->getScore(outG.getConfig())));
  const void* address = static_cast<const void*>(n);
  std::stringstream ss;
  ss << address;

  _w.writeText(ss.str());
  _w.closeTag();

}

// _____________________________________________________________________________
void SvgOutput::renderEdgeTripGeom(const graph::TransitGraph& outG,
    const graph::EdgeTripGeom& g, const graph::Edge* e, double w, double h) {

  const graph::NodeFront* nfTo = e->getTo()->getNodeFrontFor(e);
  const graph::NodeFront* nfFrom = e->getFrom()->getNodeFrontFor(e);

  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  geo::PolyLine center = g.getGeom();
  center.applyChaikinSmooth(3);

  double lineW = g.getWidth();
  double lineSpc = g.getSpacing();
  double offsetStep = lineW + lineSpc;
  double oo = g.getTotalWidth();

  double o = oo;

  assert(outG.getConfig().find(&g) != outG.getConfig().end());

  util::geo::Point dir = center.getPointAt(1).p;
  std::stringstream style;
  style.clear();
  style << "fill:none;stroke:white"
    << ";stroke-linecap:round;stroke-opacity:1;stroke-width:.5";
  printPoint(dir, style.str(), w, h, xOffset, yOffset);

  size_t a = 0;
  for (size_t i : outG.getConfig().find(&g)->second) {
      const graph::TripOccurance& r = g.getTripsUnordered()[i];
      geo::PolyLine p = center;

      double offset = -(o - oo / 2.0 - g.getWidth() /2.0);

      p.offsetPerp(offset);
      p.applyChaikinSmooth(3);
      p.simplify(0.5);

      // TODO: why is this check necessary? shouldnt be!
      // ___ OUTFACTOR
      if (nfTo && nfFrom && nfTo->geom.getLine().size() > 0 && nfFrom->geom.getLine().size() > 0) {
        if (g.getGeomDir() == e->getTo()) {
          std::set<geo::PointOnLine, geo::PointOnLineCompare> iSects = nfTo->geom.getIntersections(p);
          if (iSects.size() > 0) {
            p = p.getSegment(0, iSects.begin()->totalPos);
          } else {
            p << nfTo->geom.projectOn(p.getLine().back()).p;
          }

          std::set<geo::PointOnLine, geo::PointOnLineCompare> iSects2 = nfFrom->geom.getIntersections(p);
          if (iSects2.size() > 0) {
            p = p.getSegment(iSects2.begin()->totalPos, 1);
          } else {
            p >> nfFrom->geom.projectOn(p.getLine().front()).p;
          }
        } else {
          p << nfFrom->geom.projectOn(p.getLine().back()).p;
          p >> nfTo->geom.projectOn(p.getLine().front()).p;

          std::set<geo::PointOnLine, geo::PointOnLineCompare> iSects = nfFrom->geom.getIntersections(p);
          if (iSects.size() > 0) {
            p = p.getSegment(0, iSects.begin()->totalPos);
          } else {
            p << nfFrom->geom.projectOn(p.getLine().back()).p;
          }

          std::set<geo::PointOnLine, geo::PointOnLineCompare> iSects2 = nfTo->geom.getIntersections(p);
          if (iSects2.size() > 0) {
            p = p.getSegment(iSects2.begin()->totalPos, 1);
          } else {
            p >> nfTo->geom.projectOn(p.getLine().front()).p;
          }
        }
      }

      renderLinePart(p, lineW, *r.route);

      std::map<std::string, std::string> tparams;
      tparams["x"] = std::to_string((p.getPointAt(0.5).p.get<0>() - xOffset) * _cfg->outputResolution);
      tparams["y"] = std::to_string(h-(p.getPointAt(0.5).p.get<1>() - yOffset) * _cfg->outputResolution);
      tparams["fill"] = "white";
      tparams["font-size"] = "3px";
      _w.openTag("text", tparams);
      _w.writeText(std::to_string(a));
      _w.closeTag();

      a++;

      //break;
      o -= offsetStep;
  }
}

// _____________________________________________________________________________
void SvgOutput::renderDelegates(const graph::TransitGraph& outG,
    double w, double h) {
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  for (auto& a : _delegates) {
    _w.openTag("g");
    for (auto& pd : a.second) {
      printLine(pd.first.second, pd.first.first, w, h, xOffset, yOffset);
    }
    _w.closeTag();

    _w.openTag("g");
    for (auto& pd : a.second) {
      printLine(pd.second.second, pd.second.first, w, h, xOffset, yOffset);
    }
    _w.closeTag();
  }

}

// _____________________________________________________________________________
void SvgOutput::printPoint(const util::geo::Point& p,
													const std::string& style,
                          double w, double h, int64_t xOffs, int64_t yOffs) {
  std::map<std::string, std::string> params;
  params["cx"] = std::to_string((p.get<0>() - xOffs) * _cfg->outputResolution);
  params["cy"] = std::to_string(h-(p.get<1>() - yOffs) * _cfg->outputResolution);
  params["r"] = "2";
  params["style"] = style;
  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::printLine(const transitmapper::geo::PolyLine& l,
													const std::string& style,
                          double w, double h, int64_t xOffs, int64_t yOffs) {
	std::map<std::string, std::string> params;
  params["style"] = style;
  printLine(l, params, w, h, xOffs, yOffs);
}

// _____________________________________________________________________________
void SvgOutput::printLine(const transitmapper::geo::PolyLine& l,
													const std::map<std::string, std::string>& ps,
                          double w, double h, int64_t xOffs, int64_t yOffs) {
	std::map<std::string, std::string> params = ps;
	std::stringstream points;

	for (auto& p : l.getLine()) {
		points << " " << (p.get<0>() - xOffs) * _cfg->outputResolution << ","
      << h - (p.get<1>() - yOffs) * _cfg->outputResolution;
	}

	params["points"] = points.str();

	_w.openTag("polyline", params);
	_w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::printPolygon(const util::geo::Polygon& g,
													const std::map<std::string, std::string>& ps,
                          double w, double h, int64_t xOffs, int64_t yOffs)
{
	std::map<std::string, std::string> params = ps;
	std::stringstream points;

	for (auto& p : g.outer()) {
		points << " " << (p.get<0>() - xOffs) * _cfg->outputResolution
      << "," << h - (p.get<1>() - yOffs) * _cfg->outputResolution;
	}

	params["points"] = points.str();

	_w.openTag("polygon", params);
	_w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::printPolygon(const util::geo::Polygon& g,
													const std::string& style,
                          double w, double h, int64_t xOffs, int64_t yOffs) {
	std::map<std::string, std::string> params;
  params["style"] = style;
  printPolygon(g, params, w, h, xOffs, yOffs);
}
