[![2016 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and  orthoradial (right).](examples/render/stuttgart-example.png?raw=true)](xamples/render/stuttgart-example.png?raw=true)
*2016 Stuttgart light rail network maps generated from GTFS data, with optimal line orderings, geographically correct (left), octilinear (middle), and  orthoradial (right).*

LOOM
====

Software suite for the automated generation of geographically correct or schematic transit maps.

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

TODO