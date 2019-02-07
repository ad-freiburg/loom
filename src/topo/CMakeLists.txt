file(GLOB_RECURSE gtfs2geo_SRC *.cpp)

set(gtfs2geo_main Gtfs2GeoMain.cpp)

list(REMOVE_ITEM gtfs2geo_SRC ${gtfs2geo_main})
list(REMOVE_ITEM gtfs2geo_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(gtfs2geo ${gtfs2geo_main})
add_library(gtfs2geo_dep ${gtfs2geo_SRC})

target_include_directories(gtfs2geo_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
target_link_libraries(gtfs2geo gtfs2geo_dep util ad_cppgtfs ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GDAL_LIBRARY} -lproj -lpthread)
