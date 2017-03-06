file(GLOB_RECURSE skeletonbuilder_SRC *.cpp)

set(skeletonbuilder_main SkeletonBuilderMain.cpp)

list(REMOVE_ITEM skeletonbuilder_SRC ${skeletonbuilder_main})
list(REMOVE_ITEM skeletonbuilder_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(skeletonbuilder ${skeletonbuilder_main})
add_library(skeletonbuilder_dep ${skeletonbuilder_SRC})

target_link_libraries(skeletonbuilder skeletonbuilder_dep transitmap_dep gtfsparser ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GDAL_LIBRARY} -lproj -lpthread)
