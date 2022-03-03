include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(gtfs2graphTest TestMain.cpp)
target_link_libraries(gtfs2graphTest gtfs2graph_dep util ad_cppgtfs)
