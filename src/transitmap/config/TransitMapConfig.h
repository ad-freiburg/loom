// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
#define TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

#include <string>

namespace transitmapper {
namespace config {

struct Config {
  double lineWidth;
  double lineSpacing;

  double lineLabelSize;
  double stationLabelSize;

  std::string name;
  std::string dbgPath;

  std::string renderMethod;

  double outputResolution;
  double inputSmoothing;
  double innerGeometryPrecision;

  double outputPadding;

  double outlineWidth;
  std::string outlineColor;

  bool renderStations;
  bool renderNodeFronts;
  bool renderNodeCircles;
  bool renderNodePolygons;
  bool renderStationNames;
  bool renderEdges;
  bool renderLabels;

  bool outputStats;
  bool renderStats;
  bool renderNodeConnections;
  bool expandFronts;
  bool tightStations;

  bool renderDirMarkers;

  std::string worldFilePath;
};

}  // namespace config
}  // namespace transitmapper

#endif  // TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
