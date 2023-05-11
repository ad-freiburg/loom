// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
#define TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_

#include <string>

namespace transitmapper {
namespace config {

struct Config {
  double lineWidth = 20;
  double lineSpacing = 10;

  double lineLabelSize = 40;
  double stationLabelSize = 60;

  std::string renderMethod = "svg";

  std::string mvtPath = ".";

  bool writeStats = false;

  double outputResolution = 0.1;
  double inputSmoothing = 1;
  double innerGeometryPrecision = 3;

  double outputPadding = -1;

  double outlineWidth = 1;
  std::string outlineColor;

  bool renderStations = true;
  bool renderNodeFronts = false;
  bool renderNodeCircles = false;
  bool renderEdges = true;
  bool renderLabels = false;
  bool dontLabelDeg2 = false;
  bool fromDot = false;

  bool randomColors = false;

  bool renderNodeConnections = true;
  bool tightStations = false;

  std::vector<size_t> mvtZooms;

  bool renderDirMarkers = false;
  std::string worldFilePath;
};

}  // namespace config
}  // namespace transitmapper

#endif  // TRANSITMAP_CONFIG_TRANSITMAPCONFIG_H_
