// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename T>
void GeoJsonOutput::print(const Point<T>& p, json::Val attrs) {
  _wr.obj();
  _wr.keyVal("type", "Feature");

  _wr.key("geometry");
  _wr.obj();
  _wr.keyVal("type", "Point");
  _wr.key("coordinates");
  _wr.arr();
  _wr.val(p.getX());
  _wr.val(p.getY());
  _wr.close();
  _wr.close();
  _wr.key("properties");
  _wr.val(attrs);
  _wr.close();
}

// _____________________________________________________________________________
template <typename T>
void GeoJsonOutput::print(const Line<T>& line, json::Val attrs) {
  if (!line.size()) return;
  _wr.obj();
  _wr.keyVal("type", "Feature");

  _wr.key("geometry");
  _wr.obj();
  _wr.keyVal("type", "LineString");
  _wr.key("coordinates");
  _wr.arr();
  for (auto p : line) {
    _wr.arr();
    _wr.val(p.getX());
    _wr.val(p.getY());
    _wr.close();
  }
  _wr.close();
  _wr.close();
  _wr.key("properties");
  _wr.val(attrs);
  _wr.close();
}

// _____________________________________________________________________________
template <typename T>
void GeoJsonOutput::printLatLng(const Point<T>& p, json::Val attrs) {
  auto projP = util::geo::webMercToLatLng<double>(p.getX(), p.getY());
  print(projP, attrs);
}

// _____________________________________________________________________________
template <typename T>
void GeoJsonOutput::printLatLng(const Line<T>& line, json::Val attrs) {
  Line<T> projL;
  for (auto p : line) projL.push_back(util::geo::webMercToLatLng<double>(p.getX(), p.getY()));

  print(projL, attrs);
}
