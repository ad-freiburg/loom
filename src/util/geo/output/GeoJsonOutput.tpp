// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename T>
void GeoJsonOutput::print(const Point<T>& p, Attrs attrs) {
  _wr.obj();
  _wr.keyVal("type", "Feature");

  _wr.key("geometry");
  _wr.obj();
  _wr.keyVal("type", "Point");
  _wr.key("coordinates");
  _wr.arr();
  _wr.val(p.template get<0>());
  _wr.val(p.template get<1>());
  _wr.close();
  _wr.close();
  _wr.key("properties");
  _wr.obj(attrs);
  _wr.close();
}

// _____________________________________________________________________________
template <typename T>
void GeoJsonOutput::print(const Line<T>& line, Attrs attrs) {
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
    _wr.val(p.template get<0>());
    _wr.val(p.template get<1>());
    _wr.close();
  }
  _wr.close();
  _wr.close();
  _wr.key("properties");
  _wr.obj(attrs);
  _wr.close();
}
