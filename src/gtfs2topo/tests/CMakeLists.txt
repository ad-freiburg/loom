include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(gtfs2topoTest TestMain.cpp)
target_link_libraries(gtfs2topoTest gtfs2topo_dep util ad_cppgtfs ${Boost_LIBRARIES} -lproj)
