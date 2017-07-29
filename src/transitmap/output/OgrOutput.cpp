// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdint.h>
#include <boost/filesystem.hpp>
#include <ostream>
#include "ogrsf_frmts.h"
#include "transitmap/config/TransitMapConfig.h"
#include "transitmap/graph/Edge.h"
#include "transitmap/graph/TransitGraph.h"
#include "transitmap/output/OgrOutput.h"
#include "util/geo/PolyLine.h"
#include "util/log/Log.h"

using namespace util::geo;
using namespace transitmapper;
using namespace output;

// _____________________________________________________________________________
OgrOutput::OgrOutput(const std::string& outFolder, const config::Config* cfg)
    : _outFolder(outFolder), _cfg(cfg), _t(SHAPEFILE) {}

// _____________________________________________________________________________
void OgrOutput::print(const graph::TransitGraph& outG) {
  LOG(INFO) << "Writing raw graph to OGR file in " << _outFolder << "...";

  OGRDataSource* poDS = getDataSource(_t, _outFolder, false);

  if (poDS == 0) {
    LOG(ERROR) << "Could not init OGR writer.";
  }

  OGRLayer* edgeLayer = createOGREdgeLayer(poDS);

  if (edgeLayer == 0) {
    return;
  }

  if (_t == GEOJSON || _t == CSV) {
    // geojson doesnt support multiple layers, start new
    // datasource and write single layers

    OGRDataSource::DestroyDataSource(poDS);
    poDS = getDataSource(_t, _outFolder, true);

    if (poDS == 0) {
      LOG(ERROR) << "Could not init OGR writer.";
    }
  }

  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      addEdge(e, edgeLayer);
    }
  }

  // write the node layer
  OGRLayer* nodeLayer = createOGRNodeLayer(poDS);
  if (nodeLayer == 0) {
    return;
  }

  for (graph::Node* n : outG.getNodes()) {
    LOG(INFO) << n->getStops().size() << std::endl;
    addNode(n, nodeLayer);
  }

  OGRDataSource::DestroyDataSource(poDS);
  LOG(INFO) << "OGR written successfully\n";
}

// _____________________________________________________________________________
bool OgrOutput::addEdge(const graph::Edge* e, OGRLayer* layer) const {
  OGRFeature* edge;
  edge = OGRFeature::CreateFeature(layer->GetLayerDefn());

  edge->SetField("from",
                 boost::lexical_cast<std::string>(e->getFrom()).c_str());
  edge->SetField("to", boost::lexical_cast<std::string>(e->getTo()).c_str());

  OGRLineString geom;

  std::string wktStr = e->getGeom().getWKT();

  const char* wkt[1];
  wkt[0] = wktStr.c_str();
  geom.importFromWkt(const_cast<char**>(wkt));
  edge->SetGeometry(&geom);

  layer->CreateFeature(edge);
  OGRFeature::DestroyFeature(edge);

  // TODO: check errors
  return true;
}

// _____________________________________________________________________________
bool OgrOutput::addNode(const graph::Node* n, OGRLayer* layer) const {
  OGRFeature* node;
  node = OGRFeature::CreateFeature(layer->GetLayerDefn());

  node->SetField("id", boost::lexical_cast<std::string>(n).c_str());

  OGRPoint geom;

  std::stringstream wktStr;
  wktStr << std::setprecision(12) << boost::geometry::wkt(n->getPos());
  std::string wktString = wktStr.str();

  const char* wkt[1];
  wkt[0] = wktString.c_str();
  geom.importFromWkt(const_cast<char**>(wkt));
  node->SetGeometry(&geom);

  layer->CreateFeature(node);
  OGRFeature::DestroyFeature(node);

  // TODO: check errors
  return true;
}

// _____________________________________________________________________________
OGRLayer* OgrOutput::createOGREdgeLayer(OGRDataSource* poDS) const {
  OGRSpatialReference layerProjection;
  layerProjection.importFromProj4(_cfg->projectionString.c_str());

  // write the edge layer
  OGRLayer* edgeLayer;
  edgeLayer =
      poDS->CreateLayer("segments", &layerProjection, wkbLineString, NULL);

  if (edgeLayer == 0) {
    LOG(ERROR) << "Could not create layer 'segments'";
    return 0;
  }

  // fields always included
  OGRFieldDefn fromField("from", OFTString);
  OGRFieldDefn toField("to", OFTString);
  OGRFieldDefn routesField("routes", OFTString);

  edgeLayer->CreateField(&fromField);
  edgeLayer->CreateField(&toField);
  edgeLayer->CreateField(&routesField);

  return edgeLayer;
}

// _____________________________________________________________________________
OGRLayer* OgrOutput::createOGRNodeLayer(OGRDataSource* poDS) const {
  OGRSpatialReference layerProjection;
  layerProjection.importFromProj4(_cfg->projectionString.c_str());

  // write the edge layer
  OGRLayer* nodeLayer;
  nodeLayer = poDS->CreateLayer("nodes", &layerProjection, wkbPoint, NULL);

  if (nodeLayer == 0) {
    LOG(ERROR) << "Could not create layer 'node'";
    return 0;
  }

  // fields always included
  OGRFieldDefn idField("id", OFTString);

  nodeLayer->CreateField(&idField);

  return nodeLayer;
}

// _____________________________________________________________________________
OGRDataSource* OgrOutput::getDataSource(OGROutputType t,
                                        const std::string& folder,
                                        bool reuse) const {
  OGRSFDriver* poDriver = 0;

  OGRRegisterAll();

  switch (t) {
    case SHAPEFILE:
      poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName(
          "ESRI Shapefile");
      break;
    case GEOJSON:
      poDriver =
          OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("GeoJSON");
      break;
    case GML:
      poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("GML");
      break;
    case CSV:
      poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("CSV");
      break;
    case GMT:
      poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("GMT");
      break;
    case KML:
      poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("KML");
      break;
    case MAPINFO:
      poDriver =
          OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("MapInfo File");
      break;
    case PDF:
      poDriver = OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("PDF");
      break;
    case POSTGRES:
      poDriver =
          OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("PGDump");
      break;
    case GEOCONCEPT:
      poDriver =
          OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("Geoconcept");
      break;
    case SQLITE:
      poDriver =
          OGRSFDriverRegistrar::GetRegistrar()->GetDriverByName("SQLite");
  }

  if (poDriver == 0) {
    LOG(ERROR) << "Could not init OGR writer!";
    return 0;
  }

  boost::filesystem::path dir(folder);

  // some drivers dont create the directory, some dont like the directory
  // to already exists, catch inconsistent behavior here
  if (t == SHAPEFILE || t == GEOJSON || t == GML || t == PDF || t == POSTGRES ||
      t == GMT || t == KML || t == SQLITE || t == GEOCONCEPT) {
    try {
      boost::filesystem::create_directory(dir);
      if (!boost::filesystem::is_directory(dir)) {
        LOG(ERROR) << "could not create OGR output dir '" << dir << "'";
        return 0;
      }

      if (!reuse && (boost::filesystem::is_directory(dir) &&
                     !boost::filesystem::is_empty(dir))) {
        LOG(ERROR) << "OGR output directory '" << dir
                   << "' is not empty! Must be!";
        return 0;
      }
    } catch (boost::filesystem::filesystem_error e) {
      LOG(ERROR) << e.what();
      return 0;
    }
  }

  string ogrDir = dir.string();

  switch (t) {
    case CSV: {
      const char* opt[2];
      opt[0] = "GEOMETRY=AS_WKT";
      opt[1] = "LINEFORMAT=LF";
      OGRDataSource* ret =
          poDriver->CreateDataSource(ogrDir.c_str(), const_cast<char**>(opt));
      return ret;
    }
    default:
      return poDriver->CreateDataSource(ogrDir.c_str(), 0);
  }
}

// _____________________________________________________________________________
void OgrOutput::print(const optim::OptGraph& outG) {
  LOG(INFO) << "Writing optim graph to OGR file in " << _outFolder << "...";

  OGRDataSource* poDS = getDataSource(_t, _outFolder, false);

  if (poDS == 0) {
    LOG(ERROR) << "Could not init OGR writer.";
  }

  OGRLayer* edgeLayer = createOGREdgeLayer(poDS);

  if (edgeLayer == 0) {
    return;
  }

  if (_t == GEOJSON || _t == CSV) {
    // geojson doesnt support multiple layers, start new
    // datasource and write single layers

    OGRDataSource::DestroyDataSource(poDS);
    poDS = getDataSource(_t, _outFolder, true);

    if (poDS == 0) {
      LOG(ERROR) << "Could not init OGR writer.";
    }
  }

  for (optim::OptNode* n : outG.getNodes()) {
    for (optim::OptEdge* e : n->adjListOut) {
      addEdge(e, edgeLayer);
    }
  }

  // write the node layer
  OGRLayer* nodeLayer = createOGRNodeLayer(poDS);
  if (nodeLayer == 0) {
    return;
  }

  for (optim::OptNode* n : outG.getNodes()) {
    addNode(n, nodeLayer);
  }

  OGRDataSource::DestroyDataSource(poDS);
  LOG(INFO) << "OGR written successfully\n";
}

// _____________________________________________________________________________
bool OgrOutput::addNode(const optim::OptNode* n, OGRLayer* layer) const {
  OGRFeature* node;
  node = OGRFeature::CreateFeature(layer->GetLayerDefn());

  node->SetField("id", boost::lexical_cast<std::string>(n).c_str());

  OGRPoint geom;

  std::stringstream wktStr;
  wktStr << std::setprecision(12) << boost::geometry::wkt(n->node->getPos());
  std::string wktString = wktStr.str();

  const char* wkt[1];
  wkt[0] = wktString.c_str();
  geom.importFromWkt(const_cast<char**>(wkt));
  node->SetGeometry(&geom);

  layer->CreateFeature(node);
  OGRFeature::DestroyFeature(node);

  // TODO: check errors
  return true;
}

// _____________________________________________________________________________
bool OgrOutput::addEdge(const optim::OptEdge* e, OGRLayer* layer) const {
  OGRFeature* edge;
  edge = OGRFeature::CreateFeature(layer->GetLayerDefn());

  edge->SetField("from", boost::lexical_cast<std::string>(e->from).c_str());
  edge->SetField("to", boost::lexical_cast<std::string>(e->to).c_str());

  OGRLineString geom;

  PolyLine a;
  a << e->from->node->getPos();
  a << e->to->node->getPos();
  std::string wktStr = a.getWKT();

  const char* wkt[1];
  wkt[0] = wktStr.c_str();
  geom.importFromWkt(const_cast<char**>(wkt));
  edge->SetGeometry(&geom);

  layer->CreateFeature(edge);
  OGRFeature::DestroyFeature(edge);

  // TODO: check errors
  return true;
}
