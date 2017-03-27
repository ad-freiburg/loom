file(GLOB_RECURSE graphbuilder_SRC *.cpp)

set(graphbuilder_main GraphBuilderMain.cpp)

list(REMOVE_ITEM graphbuilder_SRC ${graphbuilder_main})
list(REMOVE_ITEM graphbuilder_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(graphbuilder ${graphbuilder_main})
add_library(graphbuilder_dep ${graphbuilder_SRC})

target_link_libraries(graphbuilder graphbuilder_dep pbutil gtfsparser ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GDAL_LIBRARY} -lproj -lpthread)
