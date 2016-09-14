// Copy//right 2016, U  //niversity of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <stdint.h>
#include <ostream>
#include "./SvgOutput.h"
#include "../geo/PolyLine.h"

using namespace transitmapper;
using namespace output;


// _____________________________________________________________________________
SvgOutput::SvgOutput(std::ostream* o, double scale)
: _o(o), _w(o, true), _scale(scale) {

}

// _____________________________________________________________________________
void SvgOutput::print(const graph::TransitGraph& outG) {
  std::map<std::string, std::string> params;

  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();


  int64_t width = outG.getBoundingBox().max_corner().get<0>() - xOffset;
  int64_t height = outG.getBoundingBox().max_corner().get<1>() - yOffset;

  width *= _scale;
  height *= _scale;

  params["width"] = std::to_string(width) + "px";
  params["height"] = std::to_string(height) + "px";

  *_o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  *_o << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">";

  _w.openTag("svg", params);

  // TODO: output edges

  outputEdges(outG, width, height);
  outputNodes(outG, width, height);


  _w.closeTags();
}

// _____________________________________________________________________________
void SvgOutput::outputNodes(const graph::TransitGraph& outG, double w, double h) {
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    renderNodeConnections(outG, n, w, h);
    renderNodeScore(outG, n, w, h);
  }
  _w.closeTag();

  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    std::map<std::string, std::string> params;

    if (n->getStops().size() > 0) {
      params["stroke"] = "black";
      params["stroke-width"] = "1";
      params["fill"] = "white";
      printPolygon(n->getConvexFrontHull(20), params, w, h, xOffset, yOffset);
    } else if (false) {
      params["r"] = "5";
      params["fill"] = "#FF00FF";
    }

    for (auto& f : n->getMainDirs()) {
      const geo::PolyLine p = f.geom;
      std::stringstream style;
      style << "fill:none;stroke:red"
        << ";stroke-linecap:round;stroke-opacity:0.5;stroke-width:1";
      std::map<std::string, std::string> params;
      params["style"] = style.str();
      //printLine(p, params, w, h, xOffset, yOffset);
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::outputEdges(const graph::TransitGraph& outG, double w, double h) {

  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (const graph::EdgeTripGeom& g : *e->getEdgeTripGeoms()) {
        renderEdgeTripGeom(outG, g, e, w, h);
      }
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::renderNodeConnections(const graph::TransitGraph& outG,
    const graph::Node* n, double w, double h) {
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  for (auto& ie : n->getInnerGeometries(true)) {
    std::stringstream style;
    style << "fill:none;stroke:#" << ie.route->getColorString()
      << ";stroke-linecap:round;stroke-opacity:1;stroke-width:" << ie.etg->getWidth() * _scale;
    std::map<std::string, std::string> params;
    params["style"] = style.str();
    printLine(ie.geom, params, w, h, xOffset, yOffset);
  }
}

// _____________________________________________________________________________
void SvgOutput::renderNodeScore(const graph::TransitGraph& outG,
    const graph::Node* n, double w, double h) {
  int64_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  int64_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  std::map<std::string, std::string> params;
  params["x"] = std::to_string((n->getPos().get<0>() - xOffset) * _scale + 10);
  params["y"] = std::to_string(h-(n->getPos().get<1>() - yOffset) * _scale - 6);
  params["style"] = "font-family:Verdana;font-size:8px; font-style:normal; font-weight: normal; fill: white; stroke-width: 0.2px; stroke-linecap: butt; stroke-linejoin: miter; stroke: black";
  _w.openTag("text", params);
  if (n->getStops().size()) {
    _w.writeText((*n->getStops().begin())->getName());
    _w.writeText("\n");
  }
  //_w.writeText(std::to_string(n->getScore()));
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
  center.simplify(1);
  double lineW = g.getWidth();
  double lineSpc = g.getSpacing();
  double oo = g.getTotalWidth();

  double o = oo;

  assert(g.getTripOrdering().size() == g.getTripsUnordered().size());

  for (size_t i : g.getTripOrdering()) {
      const graph::TripOccurance& r = g.getTripsUnordered()[i];
      geo::PolyLine p = center;
      p.offsetPerp(-(o - oo / 2 - g.getWidth() /2));

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

      // _______ /OUTFACTOR
      std::stringstream style;
      style << "fill:none;stroke:#" << r.route->getColorString()
        << ";stroke-linecap:round;stroke-opacity:1;stroke-width:" << lineW * _scale;
      std::map<std::string, std::string> params;
      std::stringstream id;
      params["style"] = style.str();
      printLine(p, params, w, h, xOffset, yOffset);

      /**
      std::map<std::string, std::string> tparams;
      tparams["x"] = std::to_string((p.getPointAt(0.5).p.get<0>() - xOffset) * _scale);
      tparams["y"] = std::to_string(h-(p.getPointAt(0.5).p.get<1>() - yOffset) * _scale);
      tparams["fill"] = "white";
      tparams["stroke"] = "red";
      _w.openTag("text", tparams);
      _w.writeText(r.route->getId());
      _w.closeTag();
      **/

      //break;
      o -= lineW + lineSpc;
  }
}

// _____________________________________________________________________________
void SvgOutput::printPoint(const util::geo::Point& p,
													const std::string& style,
                          double w, double h, int64_t xOffs, int64_t yOffs) {
  std::map<std::string, std::string> params;
  params["cx"] = std::to_string((p.get<0>() - xOffs) * _scale);
  params["cy"] = std::to_string(h-(p.get<1>() - yOffs) * _scale);
  params["r"] = "5";
  params["fill"] = "#FF00FF";
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
		points << " " << (p.get<0>() - xOffs)*_scale << "," << h - (p.get<1>() - yOffs) * _scale;
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
		points << " " << (p.get<0>() - xOffs)*_scale << "," << h - (p.get<1>() - yOffs) * _scale;
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
