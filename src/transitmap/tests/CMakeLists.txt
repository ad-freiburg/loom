include_directories(
	${TRANSITMAP_INCLUDE_DIR}
)

add_executable(transitmapTest TestMain.cpp)
target_link_libraries(transitmapTest transitmap_dep)
