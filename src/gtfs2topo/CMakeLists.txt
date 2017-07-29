file(GLOB_RECURSE gtfs2topo_SRC *.cpp)

set(gtfs2topo_main Gtfs2TopoMain.cpp)

list(REMOVE_ITEM gtfs2topo_SRC ${gtfs2topo_main})
list(REMOVE_ITEM gtfs2topo_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(gtfs2topo ${gtfs2topo_main})
add_library(gtfs2topo_dep ${gtfs2topo_SRC})

target_include_directories(gtfs2topo_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
target_link_libraries(gtfs2topo gtfs2topo_dep util ad_cppgtfs ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GDAL_LIBRARY} -lproj -lpthread)
