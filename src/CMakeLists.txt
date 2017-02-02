set(TRANSITMAP_INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(transitmap)
add_subdirectory(test)
add_subdirectory(gtfsparser)
add_subdirectory(skeletonbuilder)
