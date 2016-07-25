// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosip@informatik.uni-freiburg.de>

#include <stdint.h>
#include <ostream>
#include "./svgoutput.h"
#include "../geo/PolyLine.h"

using namespace transitmapper;
using namespace output;

const static double XSCALE = 0.5;
const static double YSCALE = 0.5;

// _____________________________________________________________________________
SvgOutput::SvgOutput(std::ostream* o) : _o(o), _w(o) {

}

// _____________________________________________________________________________
void SvgOutput::print(const graph::TransitGraph& outG) {


  std::map<std::string, std::string> params;

  uint32_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  uint32_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  uint32_t width = outG.getBoundingBox().max_corner().get<0>() - xOffset;
  uint32_t height = outG.getBoundingBox().max_corner().get<1>() - yOffset;

  width *= XSCALE;
  height *= YSCALE;

  params["width"] = std::to_string(width) + "px";
  params["height"] = std::to_string(height) + "px";

  _w.openTag("svg", params);

  // TODO: output edges

  outputEdges(outG, width, height);
  outputNodes(outG, width, height);

  _w.closeTags();
}

// _____________________________________________________________________________
void SvgOutput::outputNodes(const graph::TransitGraph& outG, double w, double h) {
  uint32_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  uint32_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    std::map<std::string, std::string> params;
    params["cx"] = std::to_string((n->getPos().get<0>() - xOffset) * XSCALE);
    params["cy"] = std::to_string(h-(n->getPos().get<1>() - yOffset) * YSCALE);
    params["r"] = "10";
    params["stroke"] = "black";
    params["stroke-width"] = "4";
    params["fill"] = "white";
    _w.openTag("circle", params);
    _w.closeTag();
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::outputEdges(const graph::TransitGraph& outG, double w, double h) {
  uint32_t xOffset = outG.getBoundingBox().min_corner().get<0>();
  uint32_t yOffset = outG.getBoundingBox().min_corner().get<1>();

  _w.openTag("g");
  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      for (const graph::EdgeTripGeom& g : e->getEdgeTripGeoms()) {
        // _________ outfactor this
        geo::PolyLine center = g.getGeom();
        double lineW = 8;
        double lineSpc = 3;
        double oo = (lineW + lineSpc) * g.getTrips().size();
        double o = oo;
        for (auto r : g.getTrips()) {
            geo::PolyLine p = center;
            p.offsetPerp(o - oo / 2);
            std::stringstream attrs;
            attrs << "fill:none;stroke:#" << r.first->getColorString()
              << ";stroke-width:" << lineW;
            printLine(p, attrs.str(), w, h, xOffset, yOffset);

            o -= lineW + lineSpc;
        }
        // __________
      }
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgOutput::printLine(const transitmapper::geo::PolyLine& l,
													const std::string& style,
                          double w, double h, uint32_t xOffs, uint32_t yOffs) {
	std::map<std::string, std::string> params;
	params["style"] = style;
	std::stringstream points;

	for (auto& p : l.getLine()) {
		points << " " << (p.get<0>() - xOffs)*XSCALE << "," << h - (p.get<1>() - yOffs) * YSCALE;
	}

	params["points"] = points.str();

	_w.openTag("polyline", params);

	_w.closeTag();
}
