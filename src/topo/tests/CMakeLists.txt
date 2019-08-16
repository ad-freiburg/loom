include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(topoTest TestMain.cpp)
target_link_libraries(topoTest topo_dep util ${Boost_LIBRARIES} -lproj)
