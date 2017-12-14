file(GLOB_RECURSE dot_SRC *.cpp)

set(dot2tg_main Dot2TgMain.cpp)

list(REMOVE_ITEM dot_SRC ${dot2tg_main})
list(REMOVE_ITEM dot_SRC TestMain.cpp)

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

add_subdirectory(tests)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(dot2tg ${dot2tg_main})
add_library(dot_dep ${dot_SRC})

target_include_directories(dot_dep PUBLIC ${PROJECT_SOURCE_DIR}/src/cppgtfs/src)
target_link_libraries(dot2tg dot_dep util transitmap_dep ${Boost_LIBRARIES} ${GLPK_LIBRARY} ${GDAL_LIBRARY} -lproj -lpthread)
