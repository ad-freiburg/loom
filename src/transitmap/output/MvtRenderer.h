// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_MVTRENDERER_H_
#define TRANSITMAP_OUTPUT_MVTRENDERER_H_

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "Renderer.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"
#include "util/xml/XmlWriter.h"

using util::Nullable;

namespace transitmapper {
namespace output {

enum VTGeomType : int {
  UNKNOWN = 0,
  POINT = 1,
  LINESTRING = 2,
  POLYGON = 3
};


struct Feature {
  util::geo::Line<double> line;
  std::string layer;
  Params params;
};

struct VTFeature {
  std::vector<uint32_t> tags;
  std::vector<uint32_t> geometry;
  VTGeomType type;
};

struct VTLayer {
  std::string name;
  std::vector<VTFeature> features;
  std::vector<std::string> keys;
  std::vector<std::string> values;
};

struct VTTile {
  std::vector<VTLayer> layers;
};

class MvtRenderer : public Renderer {
 public:
  MvtRenderer(const config::Config* cfg, size_t zoom);
  virtual ~MvtRenderer() {
    delete[] _grid;
    delete[] _grid2;
  };

  virtual void print(const shared::rendergraph::RenderGraph& outputGraph);

  void writeTiles(size_t z);

  void serializeTile(size_t x, size_t y, size_t z, VTTile* l);

  void printFeature(const util::geo::Line<double>& l, size_t z, size_t x,
                    size_t y, VTLayer* layer, Params params,
                    std::map<std::string, size_t>& keys,
                    std::map<std::string, size_t>& vals);

 private:
  const config::Config* _cfg;
  size_t _zoom;
  double _res;

  // tile grid
  uint64_t* _grid;
  uint64_t* _grid2;

  std::vector<Feature> _lineFeatures;

  std::vector<std::vector<size_t>> _lines;
  std::vector<std::pair<uint32_t, uint32_t>> _cells;

  std::vector<std::vector<size_t>> _lines2;
  std::vector<std::pair<uint32_t, uint32_t>> _cells2;

  mutable std::map<std::string, int> lineClassIds;
  mutable int lineClassId = 0;

  util::geo::Box<double> getBox(size_t z, size_t x, size_t y) const;
  uint32_t gridC(double c) const;
  void addFeature(const Feature& featuer);

  void outputNodes(const shared::rendergraph::RenderGraph& outputGraph);
  void outputEdges(const shared::rendergraph::RenderGraph& outputGraph);

  void renderEdgeTripGeom(const shared::rendergraph::RenderGraph& outG,
                          const shared::linegraph::LineEdge* e);

  void renderNodeConnections(const shared::rendergraph::RenderGraph& outG,
                             const shared::linegraph::LineNode* n);

  void renderLinePart(const util::geo::PolyLine<double> p, double width,
                      const shared::linegraph::Line& line,
                      const std::string& css, const std::string& oCss);

  void renderLinePart(const util::geo::PolyLine<double> p, double width,
                      const shared::linegraph::Line& line,
                      const std::string& css, const std::string& oCss,
                      const std::string& endMarker);

  void renderDelegates(const shared::rendergraph::RenderGraph& outG);

  void renderNodeFronts(const shared::rendergraph::RenderGraph& outG);

  std::multiset<InnerClique> getInnerCliques(
      const shared::linegraph::LineNode* n,
      std::vector<shared::rendergraph::InnerGeom> geoms, size_t level) const;

  void renderClique(const InnerClique& c,
                    const shared::linegraph::LineNode* node);

  bool isNextTo(const shared::rendergraph::InnerGeom& a,
                const shared::rendergraph::InnerGeom& b) const;
  bool hasSameOrigin(const shared::rendergraph::InnerGeom& a,
                     const shared::rendergraph::InnerGeom& b) const;

  size_t getNextPartner(const InnerClique& forGeom,
                        const std::vector<shared::rendergraph::InnerGeom>& pool,
                        size_t level) const;

  void writeTile(const std::vector<size_t>& objects, size_t cx, size_t cy,
                size_t z);

  std::string writeLayer(const VTLayer* l) const;
  std::string writeFeature(const VTFeature* l) const;

  std::string getLineClass(const std::string& id) const;

  std::string getMarkerPathMale(double w) const;
  std::string getMarkerPathFemale(double w) const;
};
}  // namespace output
}  // namespace transitmapper

#endif  // TRANSITMAP_OUTPUT_MVTRENDERER_H_
