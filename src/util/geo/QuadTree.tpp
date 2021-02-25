// _____________________________________________________________________________
template <typename V, typename T>
QuadTree<V, T>::QuadTree(size_t d, size_t c, const Box<T>& bbox)
    : _maxDepth(d), _capaFunc(c), _splFunc(_capaFunc) {
  _nds.push_back(QuadNode<T>{0, 0, bbox});
}

// _____________________________________________________________________________
template <typename V, typename T>
QuadTree<V, T>::QuadTree(size_t d, const SplitFunc<V, T>& splitF,
                         const Box<T>& bbox)
    : _maxDepth(d), _capaFunc(0), _splFunc(splitF) {
  _nds.push_back(QuadNode<T>{0, 0, bbox});
}

// _____________________________________________________________________________
template <typename V, typename T>
void QuadTree<V, T>::insert(const V& val, const Point<T>& pos) {
  if (!intersects(pos, _nds[0].bbox)) return;
  int64_t valId = _vals.size();
  _vals.push_back(QuadValue<V, T>{val, pos, -1});
  _vals[valId].nextValue = -1;
  insert(valId, 0, 0);
}

// _____________________________________________________________________________
template <typename V, typename T>
void QuadTree<V, T>::insert(int64_t vid, int64_t nid, size_t d) {
  if (!intersects(_vals[vid].point, _nds[nid].bbox)) return;

  if (d < _maxDepth && _nds[nid].numEls > -1 &&
      _splFunc(_nds[nid], _vals[vid])) {
    split(nid, d);
  }

  if (_nds[nid].numEls == -1) {
    // insert into fitting subtree
    for (size_t i = 0; i < 4; i++) insert(vid, _nds[nid].childs + i, d + 1);
  } else {
    if (_nds[nid].numEls == 0) {
      _nds[nid].childs = vid;
    } else {
      _vals[vid].nextValue = _nds[nid].childs;
      _nds[nid].childs = vid;
    }
    _nds[nid].numEls++;
  }
}

// _____________________________________________________________________________
template <typename V, typename T>
void QuadTree<V, T>::split(size_t nid, size_t d) {
  const auto& box = _nds[nid].bbox;
  T w = (box.getUpperRight().getX() - box.getLowerLeft().getX()) / T(2);

  int64_t curEl = _nds[nid].numEls > 0 ? _nds[nid].childs : -1;

  _nds[nid].numEls = -1;           // the node is now a leaf node
  _nds[nid].childs = _nds.size();  // the nodes quadrant block starts here

  // box at 0, 0
  _nds.push_back(QuadNode<T>{
      0, 0,
      Box<T>(box.getLowerLeft(), Point<T>(box.getLowerLeft().getX() + w,
                                          box.getLowerLeft().getY() + w))});
  // box at 0, 1
  _nds.push_back(QuadNode<T>{
      0, 0,
      Box<T>(Point<T>(box.getLowerLeft().getX() + w, box.getLowerLeft().getY()),
             Point<T>(box.getUpperRight().getX(),
                      box.getLowerLeft().getY() + w))});
  // box at 1,0
  _nds.push_back(QuadNode<T>{
      0, 0,
      Box<T>(Point<T>(box.getLowerLeft().getX(), box.getLowerLeft().getY() + w),
             Point<T>(box.getLowerLeft().getX() + w,
                      box.getUpperRight().getY()))});
  // box at 1,1
  _nds.push_back(QuadNode<T>{0, 0,
                             Box<T>(Point<T>(box.getLowerLeft().getX() + w,
                                             box.getLowerLeft().getY() + w),
                                    box.getUpperRight())});

  while (curEl > -1) {
    _vals[curEl].nextValue = -1;
    insert(curEl, nid, d + 1);
    curEl = _vals[curEl].nextValue;
  }
}

// _____________________________________________________________________________
template <typename V, typename T>
size_t QuadTree<V, T>::size() const {
  return _vals.size();
}

// _____________________________________________________________________________
template <typename V, typename T>
const QuadNode<T>& QuadTree<V, T>::getNd(size_t nid) const {
  return _nds[nid];
}

// _____________________________________________________________________________
template <typename V, typename T>
const std::vector<QuadNode<T>>& QuadTree<V, T>::getNds() const {
  return _nds;
}

// _____________________________________________________________________________
template <typename V, typename T>
void QuadTree<V, T>::print(std::ostream& o) const {
  util::geo::output::GeoJsonOutput out(o);
  for (const auto& nd : _nds) {
    if (nd.numEls == -1) continue; // don't print non-leaf nodes
    out.print(util::geo::convexHull(nd.bbox), json::Dict{{"elements", json::Int(nd.numEls)}});
  }
}
