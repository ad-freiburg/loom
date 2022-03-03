[![2016 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and  orthoradial (right).](examples/render/stuttgart-example-small.png?raw=true)](examples/render/stuttgart-example.png?raw=true)
*2016 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and  orthoradial (right).*

LOOM
====

Software suite for the automated generation of geographically correct or schematic transit maps.
Still alpha status with minimal documentation, so beware.


Requirements
------------

 * `cmake`
 * `gcc >= 4.9` (or `clang >= 5.0`)
 * Optional: `libglpk-dev`, `coinor-libcbc-dev`, `gurobi`


Building and Installation
-------------------------

Fetch this repository and init submodules:

```
git clone --recurse-submodules https://github.com/ad-freiburg/pfaedle
```

```
mkdir build && cd build
cmake ..
make -j
```

To (optionally) install, type
```
make install
```

Usage
=====

This suite consists of several tools:

* `gtfs2graph`, create a GeoJSON line graph from GTFS data
* `topo`, create an overlapping-free line graph from an arbitrary line graph
* `loom`, find optimal line orderings on a line graph
* `octi`, create a schematic version of a line graph
* `transitmap`, render a line graph into a map

All tools output a graph, in the GeoJSON format, to `stdout`, and expect a GeoJSON graph at `stdin`. Exceptions are `gtfs2graph`, where the input is a GTFS feed, and `transitmap`, which write SVG to `stdout`.

The `example` folder contains several overlapping-free line graphs.

To render the geographically correct Stuttgart map from above, use
```
cat examples/stuttgart.json | loom | transitmap > stuttgart.svg
```

To also render labels, use

```
cat examples/stuttgart.json | loom | transitmap -l > stuttgart-label.svg
```

To render an *octilinear* map, put the `octi` tool into the pipe:

```
cat examples/stuttgart.json | loom | octi | transitmap -l > stuttgart-octilin.svg
```

To render for example the orthoradial map from above, use a different base graph for `octi`:

```
cat examples/stuttgart.json | loom | octi -b orthoradial | transitmap -l > stuttgart-orthorad.svg
```

Line graph extraction from GTFS
-------------------------------

To extract for example a line graph for streetcars from GTFS data, use `gtfs2graph` as follows:

```
gtfs2graph -m tram freiburg.zip > freiburg.json
```

This line graph will have many overlapping edges and stations. To create an overlapping-free line graph ready for rendering, add `topo`:

```
gtfs2graph -m tram freiburg.zip | topo > freiburg.json
```

A full pipeline for creating an octilinear map of the Freiburg tram network would look like this:
```
gtfs2graph -m tram freiburg.zip | topo | loom | octi | transitmap > freiburg-tram.svg
```