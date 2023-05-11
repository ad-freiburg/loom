// Copyright 2023, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosip@informatik.uni-freiburg.de>

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
void RTree<V, G, T>::add(G<T> geom, V val) {
  Box<T> box = getBoundingBox(geom);

  T minCoords[2];
  T maxCoords[2];

  minCoords[0] = box.getLowerLeft().getX();
  minCoords[1] = box.getLowerLeft().getY();

  maxCoords[0] = box.getUpperRight().getX();
  maxCoords[1] = box.getUpperRight().getY();

  if (_valIdx.count(val)) assert(false);

  _valIdx[val] = box;

  _rtree->Insert(minCoords, maxCoords, val);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
bool RTree<V, G, T>::searchCb(V val, void* s) {
  static_cast<std::set<V>*>(s)->insert(val);
  return true;
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
void RTree<V, G, T>::get(const Box<T>& box, std::set<V>* s) const {
  T minCoords[2];
  T maxCoords[2];

  minCoords[0] = box.getLowerLeft().getX();
  minCoords[1] = box.getLowerLeft().getY();

  maxCoords[0] = box.getUpperRight().getX();
  maxCoords[1] = box.getUpperRight().getY();

  std::function<bool(const V&)> f = [s](const V& val) {
    s->insert(val);
    return true;
  };
  _rtree->Search(minCoords, maxCoords, f);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
void RTree<V, G, T>::get(const Box<T>& box, std::vector<V>* s) const {
  T minCoords[2];
  T maxCoords[2];

  minCoords[0] = box.getLowerLeft().getX();
  minCoords[1] = box.getLowerLeft().getY();

  maxCoords[0] = box.getUpperRight().getX();
  maxCoords[1] = box.getUpperRight().getY();

  std::function<bool(const V&)> f = [s](const V& val) {
    s->push_back(val);
    return true;
  };
  _rtree->Search(minCoords, maxCoords, f);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
void RTree<V, G, T>::remove(V val) {
  auto bit = _valIdx.find(val);
  if (bit == _valIdx.end()) return;

  Box<T> box = bit->second;

  T minCoords[2];
  T maxCoords[2];

  minCoords[0] = box.getLowerLeft().getX();
  minCoords[1] = box.getLowerLeft().getY();

  maxCoords[0] = box.getUpperRight().getX();
  maxCoords[1] = box.getUpperRight().getY();

  _valIdx.erase(bit);

  bool notFound = _rtree->Remove(minCoords, maxCoords, val);
  assert(!notFound);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
template <template <typename> class GG>
void RTree<V, G, T>::get(const GG<T>& geom, double d, std::set<V>* s) const {
  return get(util::geo::pad(getBoundingBox(geom), d), s);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
template <template <typename> class GG>
void RTree<V, G, T>::get(const std::vector<GG<T>>& geom, double d,
                         std::set<V>* s) const {
  return get(util::geo::pad(getBoundingBox(geom), d), s);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
template <template <typename> class GG>
void RTree<V, G, T>::get(const GG<T>& geom, double d, std::vector<V>* s) const {
  return get(util::geo::pad(getBoundingBox(geom), d), s);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
template <template <typename> class GG>
void RTree<V, G, T>::get(const std::vector<GG<T>>& geom, double d,
                         std::vector<V>* s) const {
  return get(util::geo::pad(getBoundingBox(geom), d), s);
}

// _____________________________________________________________________________
template <typename V, template <typename> class G, typename T>
void RTree<V, G, T>::getNeighbors(const V& val, double d,
                                  std::set<V>* s) const {
  auto bit = _valIdx.find(val);
  if (bit == _valIdx.end()) return;

  Box<T> box = util::geo::pad(bit->second, d);

  return get(box, s);
}
