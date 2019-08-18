// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename T>
Line<T> GeoGraphJsonOutput::createLine(const util::geo::Point<T>& a,
                                       const util::geo::Point<T>& b) {
  Line<T> ret;
  ret.push_back(a);
  ret.push_back(b);
  return ret;
}

// _____________________________________________________________________________
template <typename N, typename E>
void GeoGraphJsonOutput::print(const util::graph::Graph<N, E>& outG,
                               std::ostream& str) {
  printImpl(outG, str, false);
}

// _____________________________________________________________________________
template <typename N, typename E>
void GeoGraphJsonOutput::printLatLng(const util::graph::Graph<N, E>& outG,
                               std::ostream& str) {
  printImpl(outG, str, true);
}

// _____________________________________________________________________________
template <typename N, typename E>
void GeoGraphJsonOutput::printImpl(const util::graph::Graph<N, E>& outG,
                               std::ostream& str, bool proj) {
  GeoJsonOutput _out(str);

  // first pass, nodes
  for (util::graph::Node<N, E>* n : outG.getNds()) {
    if (!n->pl().getGeom()) continue;

    json::Dict props{{"id", util::toString(n)},
                     {"deg", util::toString(n->getDeg())},
                     {"deg_out", util::toString(n->getOutDeg())},
                     {"deg_in", util::toString(n->getInDeg())}};

    auto addProps = n->pl().getAttrs();
    props.insert(addProps.begin(), addProps.end());

    if (proj) {
      _out.printLatLng(*n->pl().getGeom(), props);
    } else {
      _out.print(*n->pl().getGeom(), props);
    }
  }

  // second pass, edges
  for (graph::Node<N, E>* n : outG.getNds()) {
    for (graph::Edge<N, E>* e : n->getAdjListOut()) {
      // to avoid double output for undirected graphs
      if (e->getFrom() != n) continue;
      json::Dict props{{"from", util::toString(e->getFrom())},
                       {"to", util::toString(e->getTo())},
                       {"id", util::toString(e)}};

      auto addProps = e->pl().getAttrs();
      props.insert(addProps.begin(), addProps.end());

      if (!e->pl().getGeom() || !e->pl().getGeom()->size()) {
        if (e->getFrom()->pl().getGeom()) {
          auto a = *e->getFrom()->pl().getGeom();
          if (e->getTo()->pl().getGeom()) {
            auto b = *e->getTo()->pl().getGeom();
            if (proj) {
              _out.printLatLng(createLine(a, b), props);
            } else {
              _out.print(createLine(a, b), props);
            }
          }
        }
      } else {
        if (proj) {
          _out.printLatLng(*e->pl().getGeom(), props);
        } else {
          _out.print(*e->pl().getGeom(), props);
        }
      }
    }
  }

  _out.flush();
}
