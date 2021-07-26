include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(loomTest TestMain.cpp)
target_link_libraries(loomTest loom_dep ${Boost_LIBRARIES})
