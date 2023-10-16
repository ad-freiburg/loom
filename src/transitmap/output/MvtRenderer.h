// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_MVTRENDERER_H_
#define TRANSITMAP_OUTPUT_MVTRENDERER_H_

#ifdef PROTOBUF_FOUND

#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "Renderer.h"
#include "shared/linegraph/Line.h"
#include "shared/rendergraph/RenderGraph.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/label/Labeller.h"
#include "transitmap/output/protobuf/vector_tile.pb.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"
#include "util/xml/XmlWriter.h"

using util::Nullable;

namespace transitmapper {
namespace output {

struct MvtLineFeature {
  util::geo::Line<double> line;
  std::string layer;
  Params params;
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

  void serializeTile(size_t x, size_t y, size_t z, vector_tile::Tile* l);

  void printFeature(const util::geo::Line<double>& l, size_t z, size_t x,
                    size_t y, vector_tile::Tile_Layer* layer, Params params,
                    std::map<std::string, size_t>& keys,
                    std::map<std::string, size_t>& vals);

 private:
  const config::Config* _cfg;
  size_t _zoom;
  double _res;

  // tile grid
  uint64_t* _grid;
  uint64_t* _grid2;

  std::vector<MvtLineFeature> _lineFeatures;

  std::vector<std::vector<size_t>> _lines;
  std::vector<std::pair<uint32_t, uint32_t>> _cells;

  std::vector<std::vector<size_t>> _lines2;
  std::vector<std::pair<uint32_t, uint32_t>> _cells2;

  mutable std::map<std::string, int> lineClassIds;
  mutable int lineClassId = 0;

  util::geo::Box<double> getBox(size_t z, size_t x, size_t y) const;
  uint32_t gridC(double c) const;
  void addFeature(const MvtLineFeature& featuer);

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

  std::string getLineClass(const std::string& id) const;

  std::string getMarkerPathMale(double w) const;
  std::string getMarkerPathFemale(double w) const;
};
}  // namespace output
}  // namespace transitmapper

#endif

#endif  // TRANSITMAP_OUTPUT_MVTRENDERER_H_
