set(TRANSITMAP_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_BINARY_DIR})

include_directories(
	${TRANSITMAP_INCLUDE_DIR}
	${GDAL_INCLUDE_DIR}
)

add_subdirectory(util)
add_subdirectory(transitmap)
add_subdirectory(cppgtfs)
add_subdirectory(gtfs2topo)
