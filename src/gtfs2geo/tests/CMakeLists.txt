include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(gtfs2geoTest TestMain.cpp)
target_link_libraries(gtfs2geoTest gtfs2geo_dep util ad_cppgtfs ${Boost_LIBRARIES} -lproj)
