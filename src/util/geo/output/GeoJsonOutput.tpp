// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
GeoJsonOutput::GeoJsonOutput() {}

// _____________________________________________________________________________
template <typename N, typename E>
void GeoJsonOutput::print(const util::graph::Graph<N, E>& outG) {
  json geoj;
  geoj["type"] = "FeatureCollection";
  geoj["features"] = json::array();

  // first pass, nodes
  for (graph::Node* n : outG.getNodes()) {
    if (!n->pl().getGeom()) continue;
    json feature;
    feature["type"] = "Feature";

    feature["geometry"]["type"] = "Point";
    std::vector<double> coords;
    coords.push_back(n->pl().getGeom()->get<0>());
    coords.push_back(n->pl().getGeom()->get<1>());
    feature["geometry"]["coordinates"] = coords;

    json::object_t props = json::object();
    props["id"] = toString(n);

    n->pl().getAttrs(props);

    feature["properties"] = props;

    geoj["features"].push_back(feature);
  }

  // second pass, edges
  for (graph::Node* n : outG.getNodes()) {
    for (graph::Edge* e : n->getAdjListOut()) {
      if (!e->pl().getGeom()) continue;
      json feature;
      feature["type"] = "Feature";
      feature["properties"]["from"] = toString(e->getFrom());
      feature["properties"]["to"] = toString(e->getTo());

      feature["geometry"]["type"] = "LineString";
      feature["geometry"]["coordinates"] = json::array();

      for (auto p : *e->pl().getGeom()) {
        std::vector<double> coords;
        coords.push_back(p.get<0>());
        coords.push_back(p.get<1>());
        feature["geometry"]["coordinates"].push_back(coords);
      }

      json::object_t props = json::object();
      props["id"] = toString(n);

      n->pl().getAttrs(props);

      feature["properties"] = props;


      geoj["features"].push_back(feature);
    }
  }

  std::cout << geoj.dump(2);
}
