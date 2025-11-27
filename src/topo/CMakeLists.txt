file(GLOB_RECURSE topo_SRC *.cpp)

set(topo_main TopoMain.cpp)

list(REMOVE_ITEM topo_SRC ${topo_main})
list(REMOVE_ITEM topo_SRC TestMain.cpp)

include_directories(
	${LOOM_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(topo ${topo_main})
add_library(topo_dep ${topo_SRC})

target_link_libraries(topo topo_dep shared_dep dot_dep pb_util_geo pb_util_graph pb_util_json pb_util)
