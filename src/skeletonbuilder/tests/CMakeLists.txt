include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(skeletonbuilderTest TestMain.cpp)
target_link_libraries(skeletonbuilderTest skeletonbuilder_dep transitmap_dep gtfsparser ${Boost_LIBRARIES} -lproj)
