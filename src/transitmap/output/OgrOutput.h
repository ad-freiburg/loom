// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_OUTPUT_OGROUTPUT_H_
#define TRANSITMAP_OUTPUT_OGROUTPUT_H_

#include <string>
#include <ostream>
#include "ogrsf_frmts.h"
#include "Output.h"
#include "pbutil/geo/PolyLine.h"
#include "pbutil/geo/Geo.h"
#include "./../config/TransitMapConfig.h"
#include "./../graph/TransitGraph.h"
#include "./../graph/Edge.h"
#include "./../optim/OptGraph.h"

namespace transitmapper {
namespace output {

class OgrOutput : public Output {


 public:
  enum OGROutputType {SHAPEFILE, GEOJSON, GML, CSV, GMT,
                      KML, MAPINFO, PDF, POSTGRES, GEOCONCEPT, SQLITE};

  OgrOutput(const std::string& outFolder, const config::Config* cfg);
  virtual ~OgrOutput() {};

  void print(const graph::TransitGraph& outG);
  void print(const optim::OptGraph& g);

 private:
  std::string _outFolder;
  const config::Config* _cfg;
  OGROutputType _t;

  OGRDataSource* getDataSource(OGROutputType t, const std::string& folder,
                               bool reuse) const;

  OGRLayer* createOGREdgeLayer(OGROutputType t, const string& folder,
      OGRDataSource* poDS) const;

  OGRLayer* createOGRNodeLayer(OGROutputType t, const string& folder,
      OGRDataSource* poDS) const;

  bool addEdge(const graph::Edge* e, OGRLayer* layer) const;
  bool addNode(const graph::Node* e, OGRLayer* layer) const;
  bool addNode(const optim::OptNode* n, OGRLayer* layer) const;
  bool addEdge(const optim::OptEdge* e, OGRLayer* layer) const;
};

}}

#endif // TRANSITMAOT_OUTPUT_OGROUTPUT_H_
