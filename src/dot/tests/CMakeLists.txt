include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(dotTest TestMain.cpp)
target_link_libraries(dotTest dot_dep util ad_cppgtfs)
