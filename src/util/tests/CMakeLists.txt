include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(utilTest TestMain.cpp)
target_link_libraries(utilTest util)
