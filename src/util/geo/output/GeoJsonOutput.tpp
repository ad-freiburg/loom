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
  for (util::graph::Node<N, E>* n : outG.getNodes()) {
    if (!n->pl().getGeom()) continue;
    json feature;
    feature["type"] = "Feature";

    feature["geometry"]["type"] = "Point";
    std::vector<double> coords;
    coords.push_back(n->pl().getGeom()->template get<0>());
    coords.push_back(n->pl().getGeom()->template get<1>());
    feature["geometry"]["coordinates"] = coords;

    json::object_t props = json::object();
    props["id"] = toString(n);

    n->pl().getAttrs(props);

    feature["properties"] = props;

    geoj["features"].push_back(feature);
  }

  // second pass, edges
  for (graph::Node<N, E>* n : outG.getNodes()) {
    for (graph::Edge<N, E>* e : n->getAdjListOut()) {
      if (!e->pl().getGeom()) {
        std::cerr << "WARNING: edge " << e << " has no geometry! Ignoring on output..." << std::endl;
        continue;
      }
      json feature;
      feature["type"] = "Feature";
      json::object_t props = json::object();
      props["from"] = toString(e->getFrom());
      props["to"] = toString(e->getTo());

      feature["geometry"]["type"] = "LineString";
      feature["geometry"]["coordinates"] = json::array();

      for (auto p : *e->pl().getGeom()) {
        std::vector<double> coords;
        coords.push_back(p.template get<0>());
        coords.push_back(p.template get<1>());
        feature["geometry"]["coordinates"].push_back(coords);
      }

      props["id"] = toString(e);

      e->pl().getAttrs(props);

      feature["properties"] = props;

      geoj["features"].push_back(feature);
    }
  }

  std::cout << geoj.dump(2);
}
