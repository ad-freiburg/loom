file(GLOB_RECURSE transitmap_SRC *.cpp)

set(transitmap_main ${CMAKE_CURRENT_SOURCE_DIR})

list(REMOVE_ITEM transitmap_SRC ${transitmap_main})

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

configure_file (
  "_config.h.in"
  "_config.h"
)

add_executable(transitmap ${transitmap_main})
add_library(transitmap_dep ${transitmap_SRC})

target_link_libraries(transitmap transitmap_dep gtfsparser ${Boost_LIBRARIES} -lproj -lpthread)
