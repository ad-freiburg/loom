// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>
#include <ostream>
#include "shared/linegraph/Route.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/output/SvgRenderer.h"
#include "transitmap/graph/RenderGraph.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"

using namespace transitmapper;
using namespace output;
using shared::linegraph::Route;
using shared::linegraph::LineNode;
using util::geo::DPoint;
using util::geo::PolyLine;
using util::geo::DPolygon;
using util::geo::Polygon;
using util::geo::LinePoint;
using util::geo::LinePointCmp;

using util::toString;

// _____________________________________________________________________________
SvgRenderer::SvgRenderer(std::ostream* o, const config::Config* cfg,
                         const optim::Scorer* scorer)
    : _o(o), _w(o, true), _cfg(cfg), _scorer(scorer) {}

// _____________________________________________________________________________
void SvgRenderer::print(const graph::RenderGraph& outG) {
  std::map<std::string, std::string> params;
  RenderParams rparams;

  double p = _cfg->outputPadding;

  auto box = util::geo::pad(outG.getBBox(), p);

  rparams.xOff = box.getLowerLeft().getX();
  rparams.yOff = box.getLowerLeft().getY();

  rparams.width = box.getUpperRight().getX() - rparams.xOff;
  rparams.height = box.getUpperRight().getY() - rparams.yOff;

  rparams.width *= _cfg->outputResolution;
  rparams.height *= _cfg->outputResolution;

  params["width"] = std::to_string(rparams.width) + "px";
  params["height"] = std::to_string(rparams.height) + "px";

  *_o << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  *_o << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
         "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">";

  if (_cfg->renderEdges) {
    outputEdges(outG, rparams);
  }
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

  for (auto n : outG.getNds()) {
    if (_cfg->renderNodeConnections) {
      renderNodeConnections(outG, n, rparams);
    }
  }

  renderDelegates(outG, rparams);

  outputNodes(outG, rparams);
  if (_cfg->renderNodeFronts) {
    renderNodeFronts(outG, rparams);
  }

  _w.closeTags();
}

// _____________________________________________________________________________
void SvgRenderer::outputNodes(const graph::RenderGraph& outG,
                              const RenderParams& rparams) {
  _w.openTag("g");
  for (auto n : outG.getNds()) {
    std::map<std::string, std::string> params;

    if (_cfg->renderStations && n->pl().getStops().size() > 0 &&
        n->pl().getMainDirs().size() > 0) {
      params["stroke"] = "black";
      params["stroke-width"] =
          util::toString((_cfg->lineWidth / 2) * _cfg->outputResolution);
      params["fill"] = "white";

      printPolygon(
          outG.getStationHull(n, (_cfg->lineSpacing + _cfg->lineWidth) * 0.8,
                              _cfg->simpleRenderForTwoEdgeNodes),
          params, rparams);
    } else if (false) {
      params["r"] = "5";
      params["fill"] = "#FF00FF";
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::renderNodeFronts(const graph::RenderGraph& outG,
                                   const RenderParams& rparams) {
  _w.openTag("g");
  for (LineNode* n : outG.getNds()) {
    std::string color = n->pl().getStops().size() > 0 ? "red" : "black";
    for (auto& f : n->pl().getMainDirs()) {
      const PolyLine<double> p = f.geom;
      std::stringstream style;
      style << "fill:none;stroke:" << color
            << ";stroke-linejoin: "
               "miter;stroke-linecap:round;stroke-opacity:0.7;stroke-width:1";
      std::map<std::string, std::string> params;
      params["style"] = style.str();
      printLine(p, params, rparams);

      DPoint a = p.getPointAt(.5).p;

      std::stringstream styleA;
      styleA << "fill:none;stroke:" << color
             << ";stroke-linejoin: "
                "miter;stroke-linecap:round;stroke-opacity:1;stroke-width:.5";
      params["style"] = styleA.str();

      printLine(PolyLine<double>(*n->pl().getGeom(), a), params, rparams);
    }
  }
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::outputEdges(const graph::RenderGraph& outG,
                              const RenderParams& rparams) {
  struct cmp {
    bool operator()(const LineNode* lhs, const LineNode* rhs) const {
      return lhs->getAdjList().size() > rhs->getAdjList().size() ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              graph::RenderGraph::getConnCardinality(lhs) >
                  graph::RenderGraph::getConnCardinality(rhs)) ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              lhs > rhs);
    }
  };

  struct cmpEdge {
    bool operator()(const shared::linegraph::LineEdge* lhs, const shared::linegraph::LineEdge* rhs) const {
      return lhs->pl().getRoutes().size() < rhs->pl().getRoutes().size() ||
             (lhs->pl().getRoutes().size() == rhs->pl().getRoutes().size() && lhs < rhs);
    }
  };

  std::set<const LineNode*, cmp> nodesOrdered;
  std::set<const shared::linegraph::LineEdge*, cmpEdge> edgesOrdered;
  for (const auto n : outG.getNds()) {
    nodesOrdered.insert(n);
  }

  std::set<const shared::linegraph::LineEdge*> rendered;

  for (const auto n : nodesOrdered) {
    edgesOrdered.insert(n->getAdjList().begin(), n->getAdjList().end());

    for (const auto* e : edgesOrdered) {
      if (rendered.insert(e).second) {
        renderEdgeTripGeom(outG, e, rparams);
      }
    }
  }
}

// _____________________________________________________________________________
void SvgRenderer::renderNodeConnections(const graph::RenderGraph& outG,
                                        const LineNode* n,
                                        const RenderParams& rparams) {
  auto geoms = outG.innerGeoms(n, outG.getConfig(),
                                       _cfg->innerGeometryPrecision);

  for (auto& clique : getInnerCliques(n, geoms, 99)) {
    renderClique(clique, n);
  }
}

// _____________________________________________________________________________
std::multiset<InnerClique> SvgRenderer::getInnerCliques(
    const shared::linegraph::LineNode* n,
    std::vector<shared::linegraph::InnerGeometry> pool, size_t level) const {
  std::multiset<InnerClique> ret;

  // start with the first geom in pool
  while (!pool.empty()) {
    InnerClique cur(n, pool.front());
    pool.erase(pool.begin());

    size_t p;
    while ((p = getNextPartner(cur, pool, level)) < pool.size()) {
      cur.geoms.push_back(pool[p]);
      pool.erase(pool.begin() + p);
    }

    ret.insert(cur);
  }

  return ret;
}

// _____________________________________________________________________________
size_t SvgRenderer::getNextPartner(
    const InnerClique& forClique,
    const std::vector<shared::linegraph::InnerGeometry>& pool,
    size_t level) const {
  for (size_t i = 0; i < pool.size(); i++) {
    const auto& ic = pool[i];
    for (auto& ciq : forClique.geoms) {
      if (isNextTo(ic, ciq) || (level > 1 && hasSameOrigin(ic, ciq))) {
        return i;
      }
    }
  }

  return pool.size();
}

// _____________________________________________________________________________
bool SvgRenderer::isNextTo(const shared::linegraph::InnerGeometry& a,
                           const shared::linegraph::InnerGeometry b) const {
  // TODO!!!!!
  return false;

  double THRESHOLD = 0.5 * M_PI + 0.1;
  if (a.from.front == b.from.front && a.to.front == b.to.front) {
    if ((a.slotFrom - b.slotFrom == 1 || b.slotFrom - a.slotFrom == 1) &&
        (a.slotTo - b.slotTo == 1 || b.slotTo - a.slotTo == 1)) {
      double ang1 = fabs(util::geo::angBetween(a.geom.front(), a.geom.back()));
      double ang2 = fabs(util::geo::angBetween(b.geom.front(), b.geom.back()));

      return ang1 > THRESHOLD && ang2 > THRESHOLD;
    }
  }

  if (a.to.front == b.from.front && a.from.front == b.to.front) {
    if ((a.slotTo - b.slotFrom == 1 || b.slotFrom - a.slotTo == 1) &&
        (a.slotFrom - b.slotTo == 1 || b.slotTo - a.slotFrom == 1)) {
      double ang1 = fabs(util::geo::angBetween(a.geom.front(), a.geom.back()));
      double ang2 = fabs(util::geo::angBetween(b.geom.front(), b.geom.back()));

      return ang1 > THRESHOLD && ang2 > THRESHOLD;
    }
  }

  return false;
}

// _____________________________________________________________________________
bool SvgRenderer::hasSameOrigin(
    const shared::linegraph::InnerGeometry& a,
    const shared::linegraph::InnerGeometry b) const {
  if (a.from.front == b.from.front) {
    return a.slotFrom == b.slotFrom;
  }
  if (a.to.front == b.from.front) {
    return a.slotTo == b.slotFrom;
  }
  if (a.to.front == b.to.front) {
    return a.slotTo == b.slotTo;
  }
  if (a.from.front == b.to.front) {
    return a.slotFrom == b.slotTo;
  }

  return false;
}

// _____________________________________________________________________________
void SvgRenderer::renderClique(const InnerClique& cc, const LineNode* n) {
  _innerDelegates.push_back(
      std::map<uintptr_t, std::vector<OutlinePrintPair>>());
  std::multiset<InnerClique> renderCliques = getInnerCliques(n, cc.geoms, 0);
  for (const auto& c : renderCliques) {
    // the longest geom will be the ref geom
    shared::linegraph::InnerGeometry ref = c.geoms[0];
    for (size_t i = 1; i < c.geoms.size(); i++) {
      if (c.geoms[i].geom.getLength() > ref.geom.getLength()) ref = c.geoms[i];
    }

    bool raw = false;
    if (ref.geom.getLength() < _cfg->lineWidth * 2) {
      raw = true;
    }

    if (ref.geom.getLength() < _cfg->lineWidth / 8) {
      // continue;
    }

    for (size_t i = 0; i < c.geoms.size(); i++) {
      PolyLine<double> pl = c.geoms[i].geom;

      if (!raw) {
        double off = -(_cfg->lineWidth + _cfg->lineSpacing) *
                     (static_cast<int>(c.geoms[i].slotFrom) -
                      static_cast<int>(ref.slotFrom));

        if (ref.from.edge->getTo() == n) off = -off;

        pl = ref.geom.offsetted(off);

        std::set<LinePoint<double>, LinePointCmp<double>> a;
        std::set<LinePoint<double>, LinePointCmp<double>> b;

        if (ref.from.front) a = ref.from.front->geom.getIntersections(pl);
        if (ref.to.front) b = ref.to.front->geom.getIntersections(pl);

        if (a.size() == 1 && b.size() == 1) {
          pl = pl.getSegment(a.begin()->totalPos, b.begin()->totalPos);
        } else if (a.size() == 1) {
          pl = pl.getSegment(a.begin()->totalPos, 1);
        } else if (b.size() == 1) {
          pl = pl.getSegment(0, b.begin()->totalPos);
        }
      }

      std::stringstream styleOutlineCropped;
      styleOutlineCropped << "fill:none;stroke:#000000";

      styleOutlineCropped << ";stroke-linecap:butt;stroke-width:"
                          << (_cfg->lineWidth + _cfg->outlineWidth) *
                                 _cfg->outputResolution;
      Params paramsOutlineCropped;
      paramsOutlineCropped["style"] = styleOutlineCropped.str();

      std::stringstream styleStr;
      styleStr << "fill:none;stroke:#" << c.geoms[i].from.route->getColor();

      styleStr << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
               << _cfg->lineWidth * _cfg->outputResolution;
      Params params;
      params["style"] = styleStr.str();

      _innerDelegates.back()[(uintptr_t)c.geoms[i].from.route].push_back(
          OutlinePrintPair(PrintDelegate(params, pl),
                           PrintDelegate(paramsOutlineCropped, pl)));
    }
  }
}

// _____________________________________________________________________________
void SvgRenderer::renderLinePart(const PolyLine<double> p, double width,
                                 const Route& route, const shared::linegraph::LineEdge* edge) {
  renderLinePart(p, width, route, edge, "");
}

// _____________________________________________________________________________
void SvgRenderer::renderLinePart(const PolyLine<double> p, double width,
                                 const Route& route, const shared::linegraph::LineEdge* edge,
                                 const std::string& endMarker) {
  if (p.getLength() < width / 2) return;
  std::stringstream styleOutline;
  styleOutline << "fill:none;stroke:#000000";

  if (!endMarker.empty()) {
    styleOutline << ";marker-end:url(#" << endMarker << "_black)";
  }

  styleOutline << ";stroke-linecap:round;stroke-width:"
               << (width + _cfg->outlineWidth) * _cfg->outputResolution;
  Params paramsOutline;
  paramsOutline["style"] = styleOutline.str();

  std::stringstream styleStr;
  styleStr << "fill:none;stroke:#" << route.getColor();

  if (!endMarker.empty()) {
    styleStr << ";marker-end:url(#" << endMarker << ")";
  }

  styleStr << ";stroke-linecap:round;stroke-opacity:1;stroke-width:"
           << width * _cfg->outputResolution;
  Params params;
  params["style"] = styleStr.str();

  _delegates[0].insert(_delegates[0].begin(),
                       OutlinePrintPair(PrintDelegate(params, p),
                                        PrintDelegate(paramsOutline, p)));
}

// _____________________________________________________________________________
void SvgRenderer::renderEdgeTripGeom(const graph::RenderGraph& outG,
                                     const shared::linegraph::LineEdge* e,
                                     const RenderParams& rparams) {
  const shared::linegraph::NodeFront* nfTo = e->getTo()->pl().getNodeFrontFor(e);
  const shared::linegraph::NodeFront* nfFrom = e->getFrom()->pl().getNodeFrontFor(e);

  PolyLine<double> center = *e->pl().getGeom();

  double lineW = _cfg->lineWidth;
  double lineSpc = _cfg->lineSpacing;
  double offsetStep = lineW + lineSpc;
  double oo = outG.getTotalWidth(e);

  double o = oo;

  assert(outG.getConfig().find(e) != outG.getConfig().end());

  size_t a = 0;
  for (size_t i : outG.getConfig().find(e)->second) {
    const auto& ro = e->pl().routeOccAtPos(i);

    const Route* route = ro.route;
    PolyLine<double> p = center;

    if (p.getLength() < 0.01) continue;

    double offset = -(o - oo / 2.0 - _cfg->lineWidth / 2.0);

    p.offsetPerp(offset);

    std::set<LinePoint<double>, LinePointCmp<double>> iSects =
        nfTo->geom.getIntersections(p);
    if (iSects.size() > 0) {
      p = p.getSegment(0, iSects.begin()->totalPos);
    } else {
      p << nfTo->geom.projectOn(p.back()).p;
    }

    std::set<LinePoint<double>, LinePointCmp<double>> iSects2 =
        nfFrom->geom.getIntersections(p);
    if (iSects2.size() > 0) {
      p = p.getSegment(iSects2.begin()->totalPos, 1);
    } else {
      p >> nfFrom->geom.projectOn(p.front()).p;
    }

    double arrowLength = (_cfg->lineWidth * 2.5);

    if (_cfg->renderDirMarkers && ro.direction != 0 &&
        center.getLength() > arrowLength * 4) {
      std::stringstream markerName;
      markerName << e << ":" << route << ":" << i;

      std::string markerPathMale = getMarkerPathMale(lineW);
      std::string markerPathFemale = getMarkerPathFemale(lineW);
      EndMarker emm(markerName.str() + "_m", "#" + route->getColor(),
                    markerPathMale, lineW, lineW);
      EndMarker emf(markerName.str() + "_f", "#" + route->getColor(),
                    markerPathFemale, lineW, lineW);

      _markers.push_back(emm);
      _markers.push_back(emf);

      PolyLine<double> firstPart =
          p.getSegmentAtDist(0, p.getLength() / 2 - arrowLength / 2);
      PolyLine<double> secondPart = p.getSegmentAtDist(
          p.getLength() / 2 + arrowLength / 2, p.getLength());

      if (ro.direction == e->getTo()) {
        renderLinePart(firstPart, lineW, *route, e, markerName.str() + "_m");
        renderLinePart(secondPart.reversed(), lineW, *route, e,
                       markerName.str() + "_f");
      } else {
        renderLinePart(firstPart, lineW, *route, e, markerName.str() + "_f");
        renderLinePart(secondPart.reversed(), lineW, *route, e,
                       markerName.str() + "_m");
      }
    } else {
      renderLinePart(p, lineW, *route, e);
    }

    a++;

    o -= offsetStep;
  }
}

// _____________________________________________________________________________
std::string SvgRenderer::getMarkerPathMale(double w) const {
  std::stringstream path;
  path << "M0,0 V1 H.5 L1.3,.5 L.5,0 Z";

  return path.str();
}

// _____________________________________________________________________________
std::string SvgRenderer::getMarkerPathFemale(double w) const {
  std::stringstream path;
  path << "M1.3,1 L.5,.5 L1.3,0 L0,0 L0,1";

  return path.str();
}

// _____________________________________________________________________________
void SvgRenderer::renderDelegates(const graph::RenderGraph& outG,
                                  const RenderParams& rparams) {
  for (auto& a : _delegates) {
    _w.openTag("g");
    for (auto& pd : a.second) {
      if (_cfg->outlineWidth > 0) {
        printLine(pd.back.second, pd.back.first, rparams);
      }
      printLine(pd.front.second, pd.front.first, rparams);
    }
    _w.closeTag();
  }

  for (auto& a : _innerDelegates) {
    _w.openTag("g");
    for (auto& b : a) {
      for (auto& pd : b.second) {
        if (_cfg->outlineWidth > 0) {
          printLine(pd.back.second, pd.back.first, rparams);
        }
        // printLine(
        //    pd.back.second.getSegmentAtDist(
        //        rparams.width * _cfg->outputResolution,
        //        pd.back.second.getLength() - rparams.width),
        //    pd.back.first, rparams);
      }
      for (auto& pd : b.second) {
        printLine(pd.front.second, pd.front.first, rparams);
      }
    }
    _w.closeTag();
  }
}

// _____________________________________________________________________________
void SvgRenderer::printPoint(const DPoint& p, const std::string& style,
                             const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["cx"] =
      std::to_string((p.getX() - rparams.xOff) * _cfg->outputResolution);
  params["cy"] = std::to_string(
      rparams.height - (p.getY() - rparams.yOff) * _cfg->outputResolution);
  params["r"] = "2";
  params["style"] = style;
  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printLine(const PolyLine<double>& l, const std::string& style,
                            const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printLine(l, params, rparams);
}

// _____________________________________________________________________________
void SvgRenderer::printLine(const PolyLine<double>& l,
                            const std::map<std::string, std::string>& ps,
                            const RenderParams& rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  for (auto& p : l.getLine()) {
    points << " " << (p.getX() - rparams.xOff) * _cfg->outputResolution << ","
           << rparams.height -
                  (p.getY() - rparams.yOff) * _cfg->outputResolution;
  }

  params["points"] = points.str();

  _w.openTag("polyline", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printPolygon(const Polygon<double>& g,
                               const std::map<std::string, std::string>& ps,
                               const RenderParams& rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  for (auto& p : g.getOuter()) {
    points << " " << (p.getX() - rparams.xOff) * _cfg->outputResolution << ","
           << rparams.height -
                  (p.getY() - rparams.yOff) * _cfg->outputResolution;
  }

  params["points"] = points.str();

  _w.openTag("polygon", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printCircle(const DPoint& center, double rad,
                              const std::string& style,
                              const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printCircle(center, rad, params, rparams);
}

// _____________________________________________________________________________
void SvgRenderer::printCircle(const DPoint& center, double rad,
                              const std::map<std::string, std::string>& ps,
                              const RenderParams& rparams) {
  std::map<std::string, std::string> params = ps;
  std::stringstream points;

  params["cx"] =
      std::to_string((center.getX() - rparams.xOff) * _cfg->outputResolution);
  params["cy"] = std::to_string(
      rparams.height - (center.getY() - rparams.yOff) * _cfg->outputResolution);
  params["r"] = std::to_string(rad * _cfg->outputResolution);

  _w.openTag("circle", params);
  _w.closeTag();
}

// _____________________________________________________________________________
void SvgRenderer::printPolygon(const DPolygon& g, const std::string& style,
                               const RenderParams& rparams) {
  std::map<std::string, std::string> params;
  params["style"] = style;
  printPolygon(g, params, rparams);
}

// _____________________________________________________________________________
size_t InnerClique::getNumBranchesIn(
    const shared::linegraph::NodeFront* front) const {
  std::set<size_t> slots;
  size_t ret = 0;
  for (const auto& ig : geoms) {
    if (ig.from.front == front && !slots.insert(ig.slotFrom).second) ret++;
    if (ig.to.front == front && !slots.insert(ig.slotTo).second) ret++;
  }

  return ret;
}

// _____________________________________________________________________________
double InnerClique::getZWeight() const {
  // more weight = more to the bottom

  double BRANCH_WEIGHT = 4;

  double ret = 0;

  ret = geoms.size();  // baseline: threads with more lines to the bottom,
                       // because they are easier to follow

  for (const auto& nf : n->pl().getMainDirs()) {
    ret -= getNumBranchesIn(&nf) * BRANCH_WEIGHT;
  }

  return ret;
}

// _____________________________________________________________________________
bool InnerClique::operator<(const InnerClique& rhs) const {
  // more weight = more to the bottom
  return getZWeight() > rhs.getZWeight();
}
