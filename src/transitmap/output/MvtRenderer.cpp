// Copyright 2023, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef PROTOBUF_FOUND

#include <stdint.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <fstream>
#include <ostream>

#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/output/MvtRenderer.h"
#include "transitmap/output/protobuf/vector_tile.pb.h"
#include "util/String.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using shared::linegraph::Line;
using shared::linegraph::LineNode;
using shared::rendergraph::InnerGeom;
using shared::rendergraph::RenderGraph;
using transitmapper::label::Labeller;
using transitmapper::output::InnerClique;
using transitmapper::output::MvtRenderer;
using util::geo::Box;
using util::geo::DPoint;
using util::geo::DPolygon;
using util::geo::LinePoint;
using util::geo::LinePointCmp;
using util::geo::Polygon;
using util::geo::PolyLine;

const static double TILE_RES = 1024;

// ground width/height of a single tile on zoom level 0
const static double WEB_MERC_EXT = 20037508.3427892;

// ground width/height is halved with each zoom level
const static int GRID_ZOOM = 14;
const static int GRID_SIZE = 1 << (GRID_ZOOM);

const static double GRID_W =
    (WEB_MERC_EXT * 2) / static_cast<double>(GRID_SIZE);

const static int GRID2_ZOOM = 7;
const static int GRID2_SIZE = 1 << (GRID2_ZOOM);

const static double GRID2_W =
    (WEB_MERC_EXT * 2) / static_cast<double>(GRID2_SIZE);

// _____________________________________________________________________________
MvtRenderer::MvtRenderer(const config::Config* cfg, size_t zoom)
    : _cfg(cfg), _zoom(zoom) {
  _grid = new size_t[GRID_SIZE * GRID_SIZE];
  for (size_t i = 0; i < GRID_SIZE * GRID_SIZE; i++) {
    _grid[i] = std::numeric_limits<uint32_t>::max();
  }

  _grid2 = new size_t[GRID2_SIZE * GRID2_SIZE];
  for (size_t i = 0; i < GRID2_SIZE * GRID2_SIZE; i++) {
    _grid2[i] = std::numeric_limits<uint32_t>::max();
  }

  _res = 156543.0 / (1 << zoom);
}

// _____________________________________________________________________________
void MvtRenderer::print(const RenderGraph& outG) {
  std::map<std::string, std::string> params;

  LOGTO(DEBUG, std::cerr) << "Rendering edges...";
  if (_cfg->renderEdges) {
    outputEdges(outG);
  }

  LOGTO(DEBUG, std::cerr) << "Rendering nodes...";
  for (auto n : outG.getNds()) {
    if (_cfg->renderNodeConnections) {
      renderNodeConnections(outG, n);
    }
  }

  LOGTO(DEBUG, std::cerr) << "Writing nodes...";
  outputNodes(outG);
  if (_cfg->renderNodeFronts) {
    renderNodeFronts(outG);
  }

  LOGTO(DEBUG, std::cerr) << "Writing tiles...";
  writeTiles(_zoom);
}

// _____________________________________________________________________________
void MvtRenderer::outputNodes(const RenderGraph& outG) {
  for (auto n : outG.getNds()) {
    std::map<std::string, std::string> params;

    if (_cfg->renderStations && n->pl().stops().size() > 0 &&
        n->pl().fronts().size() > 0) {
      params["color"] = "000";
      params["fillColor"] = "fff";
      if (n->pl().stops().size()) {
        params["stationLabel"] = n->pl().stops().front().name;
        params["stationId"] = n->pl().stops().front().id;
      }
      params["width"] = util::toString((_cfg->lineWidth / 2));

      if (n->pl().getComponent() != std::numeric_limits<uint32_t>::max())
        params["component"] = util::toString(n->pl().getComponent());

      for (const auto& geom : outG.getStopGeoms(n, _cfg->tightStations, 32)) {
        addFeature({geom.getOuter(), "stations", params});
      }
    }
  }
}

// _____________________________________________________________________________
void MvtRenderer::renderNodeFronts(const RenderGraph& outG) {
  for (auto n : outG.getNds()) {
    std::string color = n->pl().stops().size() > 0 ? "red" : "black";
    for (auto& f : n->pl().fronts()) {
      const PolyLine<double> p = f.geom;

      Params params;
      params["color"] = "ff0000";
      params["lineCap"] = "round";
      params["width"] = "2";

      addFeature({p.getLine(), "lines", params});

      DPoint a = p.getPointAt(.5).p;

      addFeature(
          {PolyLine<double>(*n->pl().getGeom(), a).getLine(), "lines", params});
    }
  }
}

// _____________________________________________________________________________
void MvtRenderer::outputEdges(const RenderGraph& outG) {
  struct cmp {
    bool operator()(const LineNode* lhs, const LineNode* rhs) const {
      return lhs->getAdjList().size() > rhs->getAdjList().size() ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              RenderGraph::getConnCardinality(lhs) >
                  RenderGraph::getConnCardinality(rhs)) ||
             (lhs->getAdjList().size() == rhs->getAdjList().size() &&
              lhs > rhs);
    }
  };

  struct cmpEdge {
    bool operator()(const shared::linegraph::LineEdge* lhs,
                    const shared::linegraph::LineEdge* rhs) const {
      return lhs->pl().getLines().size() < rhs->pl().getLines().size() ||
             (lhs->pl().getLines().size() == rhs->pl().getLines().size() &&
              lhs < rhs);
    }
  };

  std::set<const LineNode*, cmp> nodesOrdered;
  std::set<const shared::linegraph::LineEdge*, cmpEdge> edgesOrdered;
  for (auto nd : outG.getNds()) nodesOrdered.insert(nd);

  std::set<const shared::linegraph::LineEdge*> rendered;

  for (const auto n : nodesOrdered) {
    edgesOrdered.insert(n->getAdjList().begin(), n->getAdjList().end());

    for (const auto* e : edgesOrdered) {
      if (rendered.insert(e).second) renderEdgeTripGeom(outG, e);
    }
  }
}

// _____________________________________________________________________________
void MvtRenderer::renderNodeConnections(const RenderGraph& outG,
                                        const LineNode* n) {
  auto geoms = outG.innerGeoms(n, _cfg->innerGeometryPrecision * _res);

  for (auto& clique : getInnerCliques(n, geoms, 9999)) renderClique(clique, n);
}

// _____________________________________________________________________________
std::multiset<InnerClique> MvtRenderer::getInnerCliques(
    const shared::linegraph::LineNode* n, std::vector<InnerGeom> pool,
    size_t level) const {
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
size_t MvtRenderer::getNextPartner(const InnerClique& forClique,
                                   const std::vector<InnerGeom>& pool,
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
bool MvtRenderer::isNextTo(const InnerGeom& a, const InnerGeom& b) const {
  double THRESHOLD = 0.5 * M_PI + 0.1;

  if (!a.from.edge) return false;
  if (!b.from.edge) return false;
  if (!a.to.edge) return false;
  if (!b.to.edge) return false;

  auto nd = RenderGraph::sharedNode(a.from.edge, a.to.edge);

  assert(a.from.edge);
  assert(b.from.edge);
  assert(a.to.edge);
  assert(b.to.edge);

  bool aFromInv = a.from.edge->getTo() == nd;
  bool bFromInv = b.from.edge->getTo() == nd;
  bool aToInv = a.to.edge->getTo() == nd;
  bool bToInv = b.to.edge->getTo() == nd;

  int aSlotFrom = !aFromInv
                      ? a.slotFrom
                      : (a.from.edge->pl().getLines().size() - 1 - a.slotFrom);
  int aSlotTo =
      !aToInv ? a.slotTo : (a.to.edge->pl().getLines().size() - 1 - a.slotTo);
  int bSlotFrom = !bFromInv
                      ? b.slotFrom
                      : (b.from.edge->pl().getLines().size() - 1 - b.slotFrom);
  int bSlotTo =
      !bToInv ? b.slotTo : (b.to.edge->pl().getLines().size() - 1 - b.slotTo);

  if (a.from.edge == b.from.edge && a.to.edge == b.to.edge) {
    if ((aSlotFrom - bSlotFrom == 1 && bSlotTo - aSlotTo == 1) ||
        (bSlotFrom - aSlotFrom == 1 && aSlotTo - bSlotTo == 1)) {
      return true;
      double ang1 = fabs(util::geo::angBetween(a.geom.front(), a.geom.back()));
      double ang2 = fabs(util::geo::angBetween(b.geom.front(), b.geom.back()));

      return ang1 > THRESHOLD && ang2 > THRESHOLD;
    }
  }

  if (a.to.edge == b.from.edge && a.from.edge == b.to.edge) {
    if ((aSlotFrom - bSlotTo == 1 && bSlotFrom - aSlotTo == 1) ||
        (bSlotTo - aSlotFrom == 1 && aSlotTo - bSlotFrom == 1)) {
      return true;
      double ang1 = fabs(util::geo::angBetween(a.geom.front(), a.geom.back()));
      double ang2 = fabs(util::geo::angBetween(b.geom.front(), b.geom.back()));

      return ang1 > THRESHOLD && ang2 > THRESHOLD;
    }
  }

  return false;
}

// _____________________________________________________________________________
bool MvtRenderer::hasSameOrigin(const InnerGeom& a, const InnerGeom& b) const {
  if (a.from.edge == b.from.edge) {
    return a.slotFrom == b.slotFrom;
  }
  if (a.to.edge == b.from.edge) {
    return a.slotTo == b.slotFrom;
  }
  if (a.to.edge == b.to.edge) {
    return a.slotTo == b.slotTo;
  }
  if (a.from.edge == b.to.edge) {
    return a.slotFrom == b.slotTo;
  }

  return false;
}

// _____________________________________________________________________________
void MvtRenderer::renderClique(const InnerClique& cc, const LineNode* n) {
  std::multiset<InnerClique> renderCliques = getInnerCliques(n, cc.geoms, 0);
  for (const auto& c : renderCliques) {
    // the longest geom will be the ref geom
    InnerGeom ref = c.geoms[0];
    for (size_t i = 1; i < c.geoms.size(); i++) {
      if (c.geoms[i].geom.getLength() > ref.geom.getLength()) ref = c.geoms[i];
    }

    for (size_t i = 0; i < c.geoms.size(); i++) {
      PolyLine<double> pl(c.geoms[i].geom);

      if (ref.geom.getLength() >
          (_cfg->lineWidth + 2 * _cfg->outlineWidth + _cfg->lineSpacing) *
              _res * 4) {
        double off =
            -(_cfg->lineWidth + _cfg->lineSpacing + 2 * _cfg->outlineWidth) *
            _res *
            (static_cast<int>(c.geoms[i].slotFrom) -
             static_cast<int>(ref.slotFrom));

        if (ref.from.edge->getTo() == n) off = -off;

        pl = ref.geom.offsetted(off);

        if (pl.getLength() / c.geoms[i].geom.getLength() > 1.5)
          pl = c.geoms[i].geom;

        std::set<LinePoint<double>, LinePointCmp<double>> a;
        std::set<LinePoint<double>, LinePointCmp<double>> b;

        if (ref.from.edge)
          a = n->pl().frontFor(ref.from.edge)->geom.getIntersections(pl);
        if (ref.to.edge)
          b = n->pl().frontFor(ref.to.edge)->geom.getIntersections(pl);

        if (a.size() > 0 && b.size() > 0) {
          pl = pl.getSegment(a.begin()->totalPos, b.begin()->totalPos);
        } else if (a.size() > 0) {
          pl = pl.getSegment(a.begin()->totalPos, 1);
        } else if (b.size() > 0) {
          pl = pl.getSegment(0, b.begin()->totalPos);
        }
      }

      if (_cfg->outlineWidth > 0) {
        Params paramsOut;
        paramsOut["color"] = "000000";
        paramsOut["line-color"] = c.geoms[i].from.line->color();
        paramsOut["line"] = c.geoms[i].from.line->label();
        paramsOut["lineCap"] = "butt";
        paramsOut["class"] = getLineClass(c.geoms[i].from.line->id());
        paramsOut["width"] =
            util::toString((2.0 * _cfg->outlineWidth + _cfg->lineWidth));

        if (n->pl().getComponent() != std::numeric_limits<uint32_t>::max())
          paramsOut["component"] = util::toString(n->pl().getComponent());

        addFeature({pl.getLine(), "inner-connections", paramsOut});
      }

      Params params;
      params["color"] = c.geoms[i].from.line->color();
      params["line-color"] = c.geoms[i].from.line->color();
      params["line"] = c.geoms[i].from.line->label();
      params["lineCap"] = "round";
      params["class"] = getLineClass(c.geoms[i].from.line->id());
      params["width"] = util::toString(_cfg->lineWidth);

      if (n->pl().getComponent() != std::numeric_limits<uint32_t>::max())
        params["component"] = util::toString(n->pl().getComponent());

      addFeature({pl.getLine(), "inner-connections", params});
    }
  }
}

// _____________________________________________________________________________
void MvtRenderer::renderEdgeTripGeom(const RenderGraph& outG,
                                     const shared::linegraph::LineEdge* e) {
  const shared::linegraph::NodeFront* nfTo = e->getTo()->pl().frontFor(e);
  const shared::linegraph::NodeFront* nfFrom = e->getFrom()->pl().frontFor(e);

  assert(nfTo);
  assert(nfFrom);

  PolyLine<double> center(*e->pl().getGeom());

  double lineW = _cfg->lineWidth;
  double outlineW = _cfg->outlineWidth;
  double lineSpc = _cfg->lineSpacing;
  double offsetStep = (lineW + 2 * outlineW + lineSpc) * _res;
  double oo = outG.getTotalWidth(e);

  double o = oo;

  for (size_t i = 0; i < e->pl().getLines().size(); i++) {
    const auto& lo = e->pl().lineOccAtPos(i);

    const Line* line = lo.line;
    PolyLine<double> p = center;

    if (p.getLength() < 0.01) continue;

    double offset =
        -(o - oo / 2.0 - ((2 * outlineW + _cfg->lineWidth) * _res) / 2.0);

    p.offsetPerp(offset);

    auto iSects = nfTo->geom.getIntersections(p);
    if (iSects.size() > 0) {
      p = p.getSegment(0, iSects.begin()->totalPos);
    } else {
      p << nfTo->geom.projectOn(p.back()).p;
    }

    auto iSects2 = nfFrom->geom.getIntersections(p);
    if (iSects2.size() > 0) {
      p = p.getSegment(iSects2.begin()->totalPos, 1);
    } else {
      p >> nfFrom->geom.projectOn(p.front()).p;
    }

    std::string css, oCss;

    if (!lo.style.isNull()) {
      css = lo.style.get().getCss();
      oCss = lo.style.get().getOutlineCss();
    }

    if (_cfg->outlineWidth > 0) {
      Params paramsOut;
      paramsOut["color"] = "000000";
      paramsOut["line-color"] = line->color();
      paramsOut["line"] = line->label();
      paramsOut["lineCap"] = "butt";
      paramsOut["class"] = getLineClass(line->id());
      paramsOut["width"] =
          util::toString((2.0 * _cfg->outlineWidth + _cfg->lineWidth));

      if (e->pl().getComponent() != std::numeric_limits<uint32_t>::max())
        paramsOut["component"] = util::toString(e->pl().getComponent());

      addFeature({p.getLine(), "lines", paramsOut});
    }

    Params params;
    params["color"] = line->color();
    params["line-color"] = line->color();
    params["line"] = line->label();
    params["lineCap"] = "round";
    params["class"] = getLineClass(line->id());
    params["width"] = util::toString(_cfg->lineWidth);

    if (e->pl().getComponent() != std::numeric_limits<uint32_t>::max())
      params["component"] = util::toString(e->pl().getComponent());

    addFeature({p.getLine(), "lines", params});

    o -= offsetStep;
  }
}

// _____________________________________________________________________________
void MvtRenderer::addFeature(const MvtLineFeature& feature) {
  double w =
      (_cfg->lineWidth + _cfg->lineSpacing + 2 * _cfg->outlineWidth) * _res;
  const auto& box = util::geo::pad(util::geo::getBoundingBox(feature.line), w);
  size_t swX = gridC(box.getLowerLeft().getX());
  size_t swY = gridC(box.getLowerLeft().getY());

  size_t neX = gridC(box.getUpperRight().getX());
  size_t neY = gridC(box.getUpperRight().getY());

  // also check surrounding
  if (swX > 0) swX = swX - 1;
  if (swY > 0) swY = swY - 1;

  if (neX < GRID_SIZE - 1) neX = neX + 1;
  if (neY < GRID_SIZE - 2) neY = neY + 1;

  for (size_t x = swX; x <= neX && x < GRID_SIZE; x++) {
    for (size_t y = swY; y <= neY && y < GRID_SIZE; y++) {
      if (util::geo::intersects(feature.line, getBox(GRID_ZOOM, x, y))) {
        if (_grid[x * GRID_SIZE + y] == std::numeric_limits<uint32_t>::max()) {
          _grid[x * GRID_SIZE + y] = _lines.size();
          _cells.push_back({x, y});
          _lines.push_back({});
        }
        _lines[_grid[x * GRID_SIZE + y]].push_back(_lineFeatures.size());
      }
    }
  }

  swX = (WEB_MERC_EXT + box.getLowerLeft().getX()) / GRID2_W;
  swY = (WEB_MERC_EXT + box.getLowerLeft().getY()) / GRID2_W;

  neX = (WEB_MERC_EXT + box.getUpperRight().getX()) / GRID2_W;
  neY = (WEB_MERC_EXT + box.getUpperRight().getY()) / GRID2_W;

  for (size_t x = swX; x <= neX && x < GRID2_SIZE; x++) {
    for (size_t y = swY; y <= neY && y < GRID2_SIZE; y++) {
      if (util::geo::intersects(feature.line, getBox(GRID2_ZOOM, x, y))) {
        if (_grid2[x * GRID2_SIZE + y] ==
            std::numeric_limits<uint32_t>::max()) {
          _grid2[x * GRID2_SIZE + y] = _lines2.size();
          _cells2.push_back({x, y});
          _lines2.push_back({});
        }

        _lines2[_grid2[x * GRID2_SIZE + y]].push_back(_lineFeatures.size());
      }
    }
  }

  _lineFeatures.push_back(feature);
}

// _____________________________________________________________________________
uint32_t MvtRenderer::gridC(double c) const {
  return (WEB_MERC_EXT + c) / GRID_W;
}

// _____________________________________________________________________________
Box<double> MvtRenderer::getBox(size_t z, size_t x, size_t y) const {
  double gridW = ((WEB_MERC_EXT * 2) / static_cast<double>(1 << z));
  Point<double> sw(x * gridW - WEB_MERC_EXT, y * gridW - WEB_MERC_EXT);
  Point<double> ne((x + 1) * gridW - WEB_MERC_EXT,
                   (y + 1) * gridW - WEB_MERC_EXT);
  return util::geo::pad(Box<double>(sw, ne), gridW * 0.1);
}

// _____________________________________________________________________________
void MvtRenderer::printFeature(const util::geo::Line<double>& l, size_t z,
                               size_t x, size_t y,
                               vector_tile::Tile_Layer* layer, Params params,
                               std::map<std::string, size_t>& keys,
                               std::map<std::string, size_t>& vals) {
  if (l.size() < 2) return;

  // skip zero-width geometries
  if ((layer->name() == "lines" || layer->name() == "inner-connections") &&
      params.count("width") && params.find("width")->second == "0")
    return;

  double tw = (WEB_MERC_EXT * 2.0) / static_cast<double>(1 << z);

  double ox = static_cast<double>(x) * tw - WEB_MERC_EXT;
  double oy = static_cast<double>(y) * tw - WEB_MERC_EXT;

  // crop
  std::vector<util::geo::Line<double>> croppedLines;

  // pad!
  auto box = util::geo::pad(getBox(z, x, y), 50 * (tw / TILE_RES));

  if (layer->name() == "stations") {
    croppedLines.push_back(l);
  } else {
    croppedLines.push_back({});

    for (size_t i = 1; i < l.size(); i++) {
      const auto& curP = l[i];
      const auto& prevP = l[i - 1];

      if (util::geo::intersects(util::geo::LineSegment<double>{curP, prevP},
                                box)) {
        croppedLines.back().push_back(prevP);
        croppedLines.back().push_back(curP);
      } else if (croppedLines.back().size() > 0) {
        croppedLines.push_back({});
      }
    }
  }

  for (const auto& ll : croppedLines) {
    // skip point-like geometries
    if (ll.size() < 2) continue;

    auto l = util::geo::simplify(ll, 2 * tw / TILE_RES);

    // skip point-like geometries
    if (ll.size() < 2) continue;
    if (l.size() == 2 && util::geo::dist(l[0], l[1]) < tw / TILE_RES) continue;

    auto feature = layer->add_features();

    for (const auto& kv : params) {
      auto kit = keys.find(kv.first);
      auto vit = vals.find(kv.second);

      if (kit != keys.end()) {
        feature->add_tags(kit->second);
      } else {
        auto k = layer->add_keys();
        *k = kv.first;
        feature->add_tags(layer->keys_size() - 1);
        keys[kv.first] = layer->keys_size() - 1;
      }

      if (vit != vals.end()) {
        feature->add_tags(vit->second);
      } else {
        auto k = layer->add_values();
        k->set_string_value(kv.second);
        feature->add_tags(layer->values_size() - 1);
        vals[kv.second] = layer->values_size() - 1;
      }
    }

    feature->set_id(1);

    if (layer->name() == "stations") {
      feature->set_type(vector_tile::Tile_GeomType_POLYGON);
    } else {
      feature->set_type(vector_tile::Tile_GeomType_LINESTRING);
    }

    // MoveTo, 1x
    feature->add_geometry((1 & 0x7) | (1 << 3));
    int px = (l[0].getX() - ox) * (TILE_RES / tw);
    int py = TILE_RES - (l[0].getY() - oy) * (TILE_RES / tw);
    feature->add_geometry((px << 1) ^ (px >> 31));
    feature->add_geometry((py << 1) ^ (py >> 31));

    // LineTo, l.size() - 1 times
    feature->add_geometry((2 & 0x7) | ((l.size() - 1) << 3));

    for (size_t i = 1; i < l.size(); i++) {
      int dx = ((l[i].getX() - ox) * (TILE_RES / tw)) - px;
      int dy = (TILE_RES - (l[i].getY() - oy) * (TILE_RES / tw)) - py;

      px += dx;
      py += dy;

      feature->add_geometry((dx << 1) ^ (dx >> 31));
      feature->add_geometry((dy << 1) ^ (dy >> 31));
    }

    if (layer->name() == "stations") {
      // close path
      feature->add_geometry((7 & 0x7) | (1 << 3));
    }
  }
}

// _____________________________________________________________________________
std::string MvtRenderer::getLineClass(const std::string& id) const {
  auto i = lineClassIds.find(id);
  if (i != lineClassIds.end()) return "line-" + std::to_string(i->second);

  lineClassIds[id] = ++lineClassId;
  return "line-" + std::to_string(lineClassId);
}

// _____________________________________________________________________________
void MvtRenderer::serializeTile(size_t x, size_t y, size_t z,
                                vector_tile::Tile* tile) {
  std::stringstream ss;

  ss << _cfg->mvtPath << "/";

  ss << z;
  int err = mkdir(ss.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  if (err == -1) {
    if (errno != EEXIST) {
      throw std::runtime_error("Could not create output directory " + ss.str());
    }
  }

  ss << "/" << x;
  err = mkdir(ss.str().c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

  if (err == -1) {
    if (errno != EEXIST) {
      throw std::runtime_error("Could not create output directory " + ss.str());
    }
  }

  ss << "/" << ((1 << z) - 1 - y) << ".mvt";

  std::fstream fo(ss.str().c_str(),
                  std::ios::out | std::ios::trunc | std::ios::binary);

  if (!fo.is_open()) {
    throw std::runtime_error("Could not open " + ss.str());
  }

  std::string a;

  tile->SerializeToString(&a);

  fo << a;
}

// _____________________________________________________________________________
void MvtRenderer::writeTiles(size_t z) {
  if (z >= GRID_ZOOM) {
    for (auto cell : _cells) {
      auto cx = cell.first;
      auto cy = cell.second;

      for (size_t ccx = (cx << (z - GRID_ZOOM));
           ccx < ((cx + 1) << (z - GRID_ZOOM)); ccx++) {
        for (size_t ccy = (cy << (z - GRID_ZOOM));
             ccy < ((cy + 1) << (z - GRID_ZOOM)); ccy++) {
          vector_tile::Tile tile;

          auto layerInner = tile.add_layers();
          layerInner->set_version(2);
          layerInner->set_name("inner-connections");
          layerInner->set_extent(TILE_RES);
          std::map<std::string, size_t> keysInner, valsInner;

          auto layerLines = tile.add_layers();
          layerLines->set_version(2);
          layerLines->set_name("lines");
          layerLines->set_extent(TILE_RES);
          std::map<std::string, size_t> keysLines, valsLines;

          auto layerStations = tile.add_layers();
          layerStations->set_version(2);
          layerStations->set_name("stations");
          layerStations->set_extent(TILE_RES);
          std::map<std::string, size_t> keysStations, valsStations;

          for (const size_t lid : _lines[_grid[cx * GRID_SIZE + cy]]) {
            const auto& l = _lineFeatures[lid];
            if (z == GRID_ZOOM ||
                util::geo::intersects(l.line, getBox(z, ccx, ccy))) {
              if (l.layer == "lines")
                printFeature(l.line, z, ccx, ccy, layerLines, l.params,
                             keysLines, valsLines);
              if (l.layer == "inner-connections")
                printFeature(l.line, z, ccx, ccy, layerInner, l.params,
                             keysInner, valsInner);
              if (l.layer == "stations")
                printFeature(l.line, z, ccx, ccy, layerStations, l.params,
                             keysStations, valsStations);
            }
          }

          serializeTile(ccx, ccy, z, &tile);
        }
      }
    }
  } else if (z >= GRID2_ZOOM) {
    for (auto cell : _cells2) {
      auto cx = cell.first;
      auto cy = cell.second;

      for (size_t ccx = (cx << (z - GRID2_ZOOM));
           ccx < ((cx + 1) << (z - GRID2_ZOOM)); ccx++) {
        for (size_t ccy = (cy << (z - GRID2_ZOOM));
             ccy < ((cy + 1) << (z - GRID2_ZOOM)); ccy++) {
          vector_tile::Tile tile;

          auto layerInner = tile.add_layers();
          layerInner->set_version(2);
          layerInner->set_name("inner-connections");
          layerInner->set_extent(TILE_RES);
          std::map<std::string, size_t> keysInner, valsInner;

          auto layerLines = tile.add_layers();
          layerLines->set_version(2);
          layerLines->set_name("lines");
          layerLines->set_extent(TILE_RES);
          std::map<std::string, size_t> keysLines, valsLines;

          auto layerStations = tile.add_layers();
          layerStations->set_version(2);
          layerStations->set_name("stations");
          layerStations->set_extent(TILE_RES);
          std::map<std::string, size_t> keysStations, valsStations;

          for (const size_t lid : _lines2[_grid2[cx * GRID2_SIZE + cy]]) {
            const auto& l = _lineFeatures[lid];
            if (z == GRID2_ZOOM ||
                util::geo::intersects(l.line, getBox(z, ccx, ccy))) {
              if (l.layer == "lines")
                printFeature(l.line, z, ccx, ccy, layerLines, l.params,
                             keysLines, valsLines);
              if (l.layer == "inner-connections")
                printFeature(l.line, z, ccx, ccy, layerInner, l.params,
                             keysInner, valsInner);
              if (l.layer == "stations")
                printFeature(l.line, z, ccx, ccy, layerStations, l.params,
                             keysStations, valsStations);
            }
          }

          serializeTile(ccx, ccy, z, &tile);
        }
      }
    }
  } else {
    for (size_t cx = 0; cx < static_cast<size_t>(1 << z); cx++) {
      for (size_t cy = 0; cy < static_cast<size_t>(1 << z); cy++) {
        vector_tile::Tile tile;

        auto layerInner = tile.add_layers();
        layerInner->set_version(2);
        layerInner->set_name("inner-connections");
        layerInner->set_extent(TILE_RES);
        std::map<std::string, size_t> keysInner, valsInner;

        auto layerLines = tile.add_layers();
        layerLines->set_version(2);
        layerLines->set_name("lines");
        layerLines->set_extent(TILE_RES);
        std::map<std::string, size_t> keysLines, valsLines;

        auto layerStations = tile.add_layers();
        layerStations->set_version(2);
        layerStations->set_name("stations");
        layerStations->set_extent(TILE_RES);
        std::map<std::string, size_t> keysStations, valsStations;

        std::vector<size_t> objects;

        if (z == 0) {
          for (size_t i = 0; i < _lineFeatures.size(); i++) {
            objects.push_back(i);
          }
        } else {
          for (size_t ccx = (cx << (GRID2_ZOOM - z));
               ccx < ((cx + 1) << (GRID2_ZOOM - z)); ccx++) {
            for (size_t ccy = (cy << (GRID2_ZOOM - z));
                 ccy < ((cy + 1) << (GRID2_ZOOM - z)); ccy++) {
              if (_grid2[ccx * GRID2_SIZE + ccy] ==
                  std::numeric_limits<uint32_t>::max())
                continue;

              const auto& lines = _lines2[_grid2[ccx * GRID2_SIZE + ccy]];
              objects.insert(objects.end(), lines.begin(), lines.end());
            }
          }

          std::sort(objects.begin(), objects.end());
        }

        for (size_t i = 0; i < objects.size(); i++) {
          if (i > 0 && objects[i] == objects[i - 1]) continue;

          const auto& l = _lineFeatures[objects[i]];

          if (l.layer == "lines")
            printFeature(l.line, z, cx, cy, layerLines, l.params, keysLines,
                         valsLines);
          if (l.layer == "inner-connections")
            printFeature(l.line, z, cx, cy, layerInner, l.params, keysInner,
                         valsInner);
          if (l.layer == "stations")
            printFeature(l.line, z, cx, cy, layerStations, l.params,
                         keysStations, valsStations);
        }

        serializeTile(cx, cy, z, &tile);
      }
    }
  }
}

#endif
