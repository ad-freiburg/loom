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
  GeoJsonOutput _out(str);

  // first pass, nodes
  for (util::graph::Node<N, E>* n : outG.getNds()) {
    if (!n->pl().getGeom()) continue;

    json::Dict props{{"id", toString(n)},
                   {"deg", toString(n->getDeg())},
                   {"deg_out", toString(n->getOutDeg())},
                   {"deg_in", toString(n->getInDeg())}};

    auto addProps = n->pl().getAttrs();
    props.insert(addProps.begin(), addProps.end());

    _out.print(*n->pl().getGeom(), props);
  }

  // second pass, edges
  for (graph::Node<N, E>* n : outG.getNds()) {
    for (graph::Edge<N, E>* e : n->getAdjListOut()) {
      // to avoid double output for undirected graphs
      if (e->getFrom() != n) continue;
      json::Dict props{{"from", toString(e->getFrom())},
                  {"to", toString(e->getTo())},
                  {"id", toString(e)}};

      auto addProps = e->pl().getAttrs();
      props.insert(addProps.begin(), addProps.end());

      if (!e->pl().getGeom() || !e->pl().getGeom()->size()) {
        if (e->getFrom()->pl().getGeom()) {
          auto a = *e->getFrom()->pl().getGeom();
          if (e->getTo()->pl().getGeom()) {
            auto b = *e->getTo()->pl().getGeom();
            _out.print(createLine(a, b), props);
          }
        }
      } else {
        _out.print(*e->pl().getGeom(), props);
      }
    }
  }

  _out.flush();
}
