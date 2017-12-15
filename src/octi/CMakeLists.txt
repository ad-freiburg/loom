file(GLOB_RECURSE octi_SRC *.cpp)

set(octi_main OctiMain.cpp)

list(REMOVE_ITEM octi_SRC ${octi_main})
list(REMOVE_ITEM octi_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(octi ${octi_main})
add_library(octi_dep ${octi_SRC})

target_include_directories(octi_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
target_link_libraries(octi octi_dep util dot_dep transitmap_dep ad_cppgtfs ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GDAL_LIBRARY} -lproj -lpthread)
