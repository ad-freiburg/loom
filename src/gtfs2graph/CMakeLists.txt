file(GLOB_RECURSE gtfs2graph_SRC *.cpp)

set(gtfs2graph_main Gtfs2GraphMain.cpp)

list(REMOVE_ITEM gtfs2graph_SRC ${gtfs2graph_main})
list(REMOVE_ITEM gtfs2graph_SRC TestMain.cpp)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(gtfs2graph ${gtfs2graph_main})
add_library(gtfs2graph_dep ${gtfs2graph_SRC})

target_include_directories(gtfs2graph_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
target_link_libraries(gtfs2graph gtfs2graph_dep util ad_cppgtfs -lpthread)
