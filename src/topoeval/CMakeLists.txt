file(GLOB_RECURSE topoeval_SRC *.cpp)

set(topoeval_main TopoEvalMain.cpp)

list(REMOVE_ITEM topoeval_SRC ${topoeval_main})

include_directories(
	${LOOM_INCLUDE_DIR}
)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(topoeval ${topoeval_main})
add_library(topoeval_dep ${topoeval_SRC})

target_link_libraries(topoeval topoeval_dep shared_dep dot_dep util)
