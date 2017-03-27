include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(graphbuilderTest TestMain.cpp)
target_link_libraries(graphbuilderTest graphbuilder_dep pbutil gtfsparser ${Boost_LIBRARIES} -lproj)
